import pysam
from cyvcf2 import VCF, Writer
import pandas as pd
import os
import vafator
import datetime
import json


class Annotator(object):

    vafator_header = {
        "name": "vafator",
        "version": vafator.VERSION,
        "date": datetime.datetime.now().ctime(),
        "timestamp": datetime.datetime.now().timestamp(),
    }

    def __init__(self, input_vcf, output_vcf, normal_bams=[], tumor_bams=[],
                 mapping_qual_thr=0, base_call_qual_thr=29):
        """
        :param input_vcf: the input VCF file
        :param normal_bams: none or many normal BAM files
        :param tumor_bams: none or many tumor BAM files
        :param output_vcf: the file path of the output VCF
        :param mapping_qual_thr: reads with a mapping quality lower or equal than this threshold will be filtered out
        :param base_call_qual_thr: bases with a basecall quality lower than or equal this threshold will be filtered out
        """
        self.mapping_quality_threshold = mapping_qual_thr
        self.base_call_quality_threshold = base_call_qual_thr
        self.vcf = VCF(input_vcf)
        # sets a line in the header with the command used to annotate the file
        self.vafator_header["input_vcf"] = input_vcf
        self.vafator_header["output_vcf"] = output_vcf
        self.vafator_header["normal_bams"] = normal_bams
        self.vafator_header["tumor_bams"] = tumor_bams
        self.vafator_header["mapping_quality_threshold"] = mapping_qual_thr
        self.vafator_header["base_call_quality_threshold"] = base_call_qual_thr
        self.vcf.add_to_header("##vafator_command_line={}".format(json.dumps(self.vafator_header)))
        # adds to the header all the names of the annotations
        for a in Annotator._get_headers("tumor", tumor_bams) + Annotator._get_headers("normal", normal_bams):
            self.vcf.add_info_to_header(a)
        self.vcf_writer = Writer(output_vcf, self.vcf)
        self.tumor_bams = [pysam.AlignmentFile(b, "rb") for b in tumor_bams]
        self.normal_bams = [pysam.AlignmentFile(b, "rb") for b in normal_bams]

    @staticmethod
    def _get_headers(prefix, bams):
        headers = []
        if len(bams) > 0:
            headers = [
                {'ID': "{}_af".format(prefix),
                 'Description': "Allele frequency for the alternate alleles in the {} samples".format(prefix),
                 'Type': 'Float', 'Number': 'A'},
                {'ID': "{}_dp".format(prefix),
                 'Description': "Total depth of coverage in the {} samples".format(prefix),
                 'Type': 'Float', 'Number': 'A'},
                {'ID': "{}_ac".format(prefix),
                 'Description': "Allele count for the alternate alleles in the {} samples".format(prefix),
                 'Type': 'Integer', 'Number': 'A'}
            ]
        if len(bams) > 1:
            for i, b in enumerate(bams, start=1):
                n = os.path.basename(b).split(".")[0]
                headers = headers + [
                    {'ID': "{}_af_{}".format(prefix, i),
                     'Description': "Allele frequency for the alternate alleles in the {} sample {}".format(prefix, n),
                     'Type': 'Float', 'Number': 'A'},
                    {'ID': "{}_dp_{}".format(prefix, i),
                     'Description': "Depth of coverage in the {} sample {}".format(prefix, n),
                     'Type': 'Float', 'Number': 'A'},
                    {'ID': "{}_ac_{}".format(prefix, i),
                     'Description': "Allele count for the alternate alleles in the {} sample {}".format(prefix, n),
                     'Type': 'Integer', 'Number': 'A'}
                ]
        return headers

    @staticmethod
    def _initialize_empty_count(bases_counts, base):
        if bases_counts.get(base) is None:
            bases_counts[base] = 0
        return bases_counts

    def _summarize_pileup(self, pileup):
        """
        :param pileup: the pileup column at a given position
        :type: pd.DataFrame
        :return: the counts in the pileup for each base including N
        :rtype: pd.Series
        """
        # filters reads by mapping quality and base call quality
        pileup = pileup[pileup.base_call_qualities > self.base_call_quality_threshold]
        pileup = pileup[pileup.mapping_qualities > self.mapping_quality_threshold]
        # ignore case in the bases as BAM specification states that it has no meaning
        bases_counts = pileup.bases.transform(lambda x: x.upper()).value_counts()
        bases_counts = Annotator._initialize_empty_count(bases_counts, 'A')
        bases_counts = Annotator._initialize_empty_count(bases_counts, 'C')
        bases_counts = Annotator._initialize_empty_count(bases_counts, 'G')
        bases_counts = Annotator._initialize_empty_count(bases_counts, 'T')
        bases_counts = Annotator._initialize_empty_count(bases_counts, 'N')
        return bases_counts

    def _get_pileup(self, variant, bam):
        """
        :param variant: the variant that determines the genomic position
        :type: cyvcf2.cyvcf2.Variant
        :param bam: the BAM file from where to read the pileups
        :type: pysam.AlignmentFile
        :return: the pileup columns at the variant position
        :rtype: pd.DataFrame
        """
        chromosome = variant.CHROM
        position = variant.POS
        # this function returns the pileups at all positions covered by reads covered the queried position
        # approximately +- read size bp
        pileups = bam.pileup(contig=chromosome, start=position-1, stop=position)
        pileup_position = -1
        while pileup_position + 1 != position:
            try:
                pileup = next(pileups)
            except StopIteration:
                return pd.DataFrame({
                    'bases': [],
                    'base_call_qualities': [],
                    'mapping_qualities': []
                })
            pileup_position = pileup.pos
        return pd.DataFrame({
            'bases': pileup.get_query_sequences(),
            'base_call_qualities': pileup.get_query_qualities(),
            'mapping_qualities': pileup.get_mapping_qualities()
        })

    @staticmethod
    def _get_af(base_counts, variant):
        """
        :type: pd.Series
        :type: cyvcf2.cyvcf2.Variant
        :return: a list containing the allele frequencies for each of the alternate alleles in the variant
        :rtype: list
        """
        return [str(base_counts.get(a, 0) / int(base_counts.sum()))
                if base_counts.sum() > 0 else str(0.0) for a in variant.ALT]

    @staticmethod
    def _get_ac(base_counts, variant):
        """
        :type: pd.Series
        :type: cyvcf2.cyvcf2.Variant
        :return: a list containing the allele counts for each of the alternate alleles in the variant
        :rtype: list
        """
        return [str(base_counts.get(a, 0)) for a in variant.ALT]

    def _write_batch(self, batch):
        for v in batch:
            self.vcf_writer.write_record(v)

    def _add_stats(self, bams,  variant, prefix):
        base_counts = [self._summarize_pileup(self._get_pileup(variant, b)) for b in bams]
        aggregated_base_counts = sum(base_counts)
        variant.INFO["{}_af".format(prefix)] = ",".join(Annotator._get_af(aggregated_base_counts, variant))
        variant.INFO["{}_ac".format(prefix)] = ",".join(Annotator._get_ac(aggregated_base_counts, variant))
        variant.INFO["{}_dp".format(prefix)] = str(aggregated_base_counts.sum())
        if len(base_counts) > 1:
            for i, c in enumerate(base_counts, start=1):
                variant.INFO["{}_af_{}".format(prefix, i)] = ",".join(Annotator._get_af(c, variant))
                variant.INFO["{}_ac_{}".format(prefix, i)] = ",".join(Annotator._get_ac(c, variant))
                variant.INFO["{}_dp_{}".format(prefix, i)] = str(c.sum())

    def run(self):
        batch = []
        for variant in self.vcf:
            # gets the counts of all bases across all BAMs
            if self.tumor_bams:
                self._add_stats(self.tumor_bams, variant, "tumor")
            if self.normal_bams:
                self._add_stats(self.normal_bams, variant, "normal")
            batch.append(variant)
            if len(batch) >= 1000:
                self._write_batch(batch)
                batch = []
        if len(batch) > 0:
            self._write_batch(batch)

        self.vcf_writer.close()
        self.vcf.close()
        for b in self.tumor_bams + self.normal_bams:
            b.close()
