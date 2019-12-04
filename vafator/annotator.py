import pysam
from cyvcf2 import VCF, Writer
import pandas as pd


class Annotator(object):

    annotations_header = [
        {'ID': 'normal_af', 'Description': 'Allele frequency for the alternate alleles in the normal sample',
         'Type': 'Float', 'Number': 'A'},
        {'ID': 'tumor_af', 'Description': 'Allele frequency for the alternate alleles in the tumor sample',
         'Type': 'Float', 'Number': 'A'},
        {'ID': 'normal_ac', 'Description': 'Allele count for the alternate alleles in the normal sample',
         'Type': 'Integer', 'Number': 'A'},
        {'ID': 'tumor_ac', 'Description': 'Allele count for the alternate alleles in the tumor sample',
         'Type': 'Integer', 'Number': 'A'},
        {'ID': 'normal_dp', 'Description': 'Depth of coverage in the normal sample', 'Type': 'Integer',
         'Number': '1'},
        {'ID': 'tumor_dp', 'Description': 'Depth of coverage in the tumor sample', 'Type': 'Integer',
         'Number': '1'}
    ]

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
        for a in self.annotations_header:
            self.vcf.add_info_to_header(a)
        self.vcf_writer = Writer(output_vcf, self.vcf)
        self.tumor_bams = [pysam.AlignmentFile(b, "rb" ) for b in tumor_bams]
        self.normal_bams = [pysam.AlignmentFile(b, "rb") for b in normal_bams]

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

    def run(self):
        for variant in self.vcf:
            # gets the counts of all bases across all BAMs
            if self.tumor_bams:
                tumor_base_counts = sum([self._summarize_pileup(self._get_pileup(variant, b)) for b in self.tumor_bams])
                variant.INFO["tumor_af"] = ",".join(Annotator._get_af(tumor_base_counts, variant))
                variant.INFO["tumor_ac"] = ",".join(Annotator._get_ac(tumor_base_counts, variant))
                variant.INFO["tumor_dp"] = str(tumor_base_counts.sum())
            if self.normal_bams:
                normal_base_counts = sum([self._summarize_pileup(self._get_pileup(variant, b)) for b in self.normal_bams])
                variant.INFO["normal_af"] = ",".join(Annotator._get_af(normal_base_counts, variant))
                variant.INFO["normal_ac"] = ",".join(Annotator._get_ac(normal_base_counts, variant))
                variant.INFO["normal_dp"] = str(normal_base_counts.sum())
            self.vcf_writer.write_record(variant)

        self.vcf_writer.close()
        self.vcf.close()