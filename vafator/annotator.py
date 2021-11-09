from typing import List
import pysam
from cyvcf2 import VCF, Writer, Variant
import pandas as pd
import os
from pandas import DataFrame, Series
from pysam import AlignmentFile
from pysam.libcalignmentfile import IteratorColumnRegion
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

    def __init__(self, input_vcf: str, output_vcf: str, normal_bams=[], tumor_bams=[],
                 mapping_qual_thr=0, base_call_qual_thr=29, prefix=None):

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
        self.vafator_header["prefix"] = prefix
        self.vcf.add_to_header("##vafator_command_line={}".format(json.dumps(self.vafator_header)))
        # adds to the header all the names of the annotations
        self.tumor_prefix = "{prefix}_tumor".format(prefix=prefix) if prefix is not None else "tumor"
        self.normal_prefix = "{prefix}_normal".format(prefix=prefix) if prefix is not None else "normal"
        for a in Annotator._get_headers(self.tumor_prefix, tumor_bams) + \
                 Annotator._get_headers(self.normal_prefix, normal_bams):
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
                 'Description': "Total depth of coverage in the {} samples (independent of alleles)".format(prefix),
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
                     'Description': "Depth of coverage in the {} sample {} (independent of alleles)".format(prefix, n),
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

    def _summarize_pileup(self, pileup: DataFrame) -> Series:

        # ignore case in the bases as BAM specification states that it has no meaning
        bases_counts = pileup.bases.transform(lambda x: x.upper()).value_counts()
        bases_counts = Annotator._initialize_empty_count(bases_counts, 'A')
        bases_counts = Annotator._initialize_empty_count(bases_counts, 'C')
        bases_counts = Annotator._initialize_empty_count(bases_counts, 'G')
        bases_counts = Annotator._initialize_empty_count(bases_counts, 'T')
        bases_counts = Annotator._initialize_empty_count(bases_counts, 'N')
        return bases_counts

    def _parse_pileup(self, pileups: IteratorColumnRegion):
        try:
            pileup = next(pileups)
        except StopIteration:
            pileup = None
        if pileup is None:
            return pd.DataFrame({
                'bases': [],
                'base_call_qualities': [],
                'mapping_qualities': []
            })
        return pd.DataFrame({
            'bases': pileup.get_query_sequences(),
            'base_call_qualities': pileup.get_query_qualities(),
            'mapping_qualities': pileup.get_mapping_qualities()
        })

    def _get_variant_pileup(self, variant: Variant, bam: AlignmentFile) -> IteratorColumnRegion:

        chromosome = variant.CHROM
        position = variant.POS
        # this function returns the pileups at all positions covered by reads covered the queried position
        # approximately +- read size bp
        return bam.pileup(contig=chromosome, start=position-1, stop=position, truncate=True,
                          min_base_quality=self.base_call_quality_threshold,
                          min_mapping_quality=self.mapping_quality_threshold)


    @staticmethod
    def _get_af(base_counts: Series, variant: Variant) -> List:

        return [str(base_counts.get(a, 0) / int(base_counts.sum()))
                if base_counts.sum() > 0 else str(0.0) for a in variant.ALT]

    @staticmethod
    def _get_ac(base_counts: Series, variant: Variant) -> List:

        return [str(base_counts.get(a, 0)) for a in variant.ALT]

    def _write_batch(self, batch):
        for v in batch:
            self.vcf_writer.write_record(v)

    def _add_snv_stats(self, bams: List[AlignmentFile], variant: Variant, prefix: str):

        base_counts = [self._summarize_pileup(self._parse_pileup(self._get_variant_pileup(variant, b))) for b in bams]
        aggregated_base_counts = sum(base_counts)
        variant.INFO["{}_af".format(prefix)] = ",".join(Annotator._get_af(aggregated_base_counts, variant))
        variant.INFO["{}_ac".format(prefix)] = ",".join(Annotator._get_ac(aggregated_base_counts, variant))
        variant.INFO["{}_dp".format(prefix)] = str(aggregated_base_counts.sum())
        if len(base_counts) > 1:
            for i, c in enumerate(base_counts, start=1):
                variant.INFO["{}_af_{}".format(prefix, i)] = ",".join(Annotator._get_af(c, variant))
                variant.INFO["{}_ac_{}".format(prefix, i)] = ",".join(Annotator._get_ac(c, variant))
                variant.INFO["{}_dp_{}".format(prefix, i)] = str(c.sum())

    def _add_insertion_stats(self, bams: List[AlignmentFile], variant: Variant, prefix: str):

        position = variant.POS
        insertion_length = len(variant.ALT[0]) - len(variant.REF)
        global_dp = 0
        global_ac = 0
        for i, b in enumerate(bams):
            pileups = self._get_variant_pileup(variant=variant, bam=b)
            try:
                pileup = next(pileups)
                dp = 0
                ac = 0
                for r in pileup.pileups:
                    dp += 1
                    if r.indel > 0:
                        # read with an insertion
                        start = r.alignment.reference_start
                        for cigar_type, cigar_length in r.alignment.cigartuples:
                            if cigar_type in [0, 2, 3, 7, 8]:  # consumes reference M, D, N, =, X
                                start += cigar_length
                                if start > position:
                                    break
                            elif cigar_type == 1:    # does not count I
                                if start == position and cigar_length == insertion_length:
                                    ac += 1
                af = self._calculate_af(ac, dp)
                if len(bams) > 1:
                    variant.INFO["{}_af_{}".format(prefix, i+1)] = af
                    variant.INFO["{}_ac_{}".format(prefix, i+1)] = ac
                    variant.INFO["{}_dp_{}".format(prefix, i+1)] = dp
                global_ac += ac
                global_dp += dp
            except StopIteration:
                # no reads
                pass

        global_af = self._calculate_af(global_ac, global_dp)
        variant.INFO["{}_af".format(prefix)] = global_af
        variant.INFO["{}_ac".format(prefix)] = global_ac
        variant.INFO["{}_dp".format(prefix)] = global_dp

    def _calculate_af(self, ac, dp):
        return float(ac) / dp if dp > 0 else 0.0

    def _add_deletion_stats(self, bams: List[AlignmentFile], variant: Variant, prefix: str):

        position = variant.POS
        deletion_length = len(variant.REF) - len(variant.ALT[0])
        global_dp = 0
        global_ac = 0
        for i, b in enumerate(bams):
            pileups = self._get_variant_pileup(variant=variant, bam=b)
            try:
                pileup = next(pileups)
                dp = 0
                ac = 0
                for r in pileup.pileups:
                    dp += 1
                    if r.indel < 0:
                        # read with an deletion
                        start = r.alignment.reference_start
                        for cigar_type, cigar_length in r.alignment.cigartuples:
                            if cigar_type in [0, 3, 7, 8]:  # consumes reference M, N, =, X
                                start += cigar_length
                            elif cigar_type == 2:   # D
                                if start == position and cigar_length == deletion_length:
                                    ac += 1
                                else:
                                    start += cigar_length
                            if start > position:
                                break
                af = self._calculate_af(ac, dp)
                if len(bams) > 1:
                    variant.INFO["{}_af_{}".format(prefix, i+1)] = af
                    variant.INFO["{}_ac_{}".format(prefix, i+1)] = ac
                    variant.INFO["{}_dp_{}".format(prefix, i+1)] = dp
                global_ac += ac
                global_dp += dp
            except StopIteration:
                # no reads
                pass

        global_af = self._calculate_af(global_ac, global_dp)
        variant.INFO["{}_af".format(prefix)] = global_af
        variant.INFO["{}_ac".format(prefix)] = global_ac
        variant.INFO["{}_dp".format(prefix)] = global_dp

    def run(self):
        batch = []
        variant: Variant
        for variant in self.vcf:
            # gets the counts of all bases across all BAMs
            if self.tumor_bams:
                if variant.is_snp:
                    self._add_snv_stats(self.tumor_bams, variant, self.tumor_prefix)
                elif variant.is_indel and not variant.is_deletion:
                    self._add_insertion_stats(self.tumor_bams, variant, self.tumor_prefix)
                elif variant.is_indel and variant.is_deletion:
                    self._add_deletion_stats(self.tumor_bams, variant, self.tumor_prefix)
            if self.normal_bams:
                if variant.is_snp:
                    self._add_snv_stats(self.normal_bams, variant, self.normal_prefix)
                elif variant.is_indel and not variant.is_deletion:
                    self._add_insertion_stats(self.normal_bams, variant, self.normal_prefix)
                elif variant.is_indel and variant.is_deletion:
                    self._add_deletion_stats(self.normal_bams, variant, self.normal_prefix)
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
