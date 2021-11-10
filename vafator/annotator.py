from typing import List
import pysam
from cyvcf2 import VCF, Writer, Variant
import os
from pysam import AlignmentFile
from pysam.libcalignmentfile import IteratorColumnRegion
import vafator
import datetime
import json

from vafator.pileups import get_insertion_metrics, get_deletion_metrics, get_snv_metrics


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

    def _get_variant_pileup(self, variant: Variant, bam: AlignmentFile) -> IteratorColumnRegion:

        chromosome = variant.CHROM
        position = variant.POS
        # this function returns the pileups at all positions covered by reads covered the queried position
        # approximately +- read size bp
        return bam.pileup(contig=chromosome, start=position-1, stop=position, truncate=True,
                          min_base_quality=self.base_call_quality_threshold,
                          min_mapping_quality=self.mapping_quality_threshold)

    def _write_batch(self, batch):
        for v in batch:
            self.vcf_writer.write_record(v)

    def _add_snv_stats(self, bams: List[AlignmentFile], variant: Variant, prefix: str):

        global_dp = 0
        global_ac = {}
        for i, b in enumerate(bams):
            pileups = self._get_variant_pileup(variant=variant, bam=b)
            ac, dp = get_snv_metrics(variant=variant, pileups=pileups)
            if len(bams) > 1:
                variant.INFO["{}_af_{}".format(prefix, i + 1)] = ",".join(
                    [str(self._calculate_af(ac[alt], dp)) for alt in variant.ALT])
                variant.INFO["{}_ac_{}".format(prefix, i + 1)] = ",".join([str(ac[alt]) for alt in variant.ALT])
                variant.INFO["{}_dp_{}".format(prefix, i + 1)] = dp
            for alt in variant.ALT:
                global_ac[alt] = global_ac.get(alt, 0) + ac[alt]
            global_dp += dp

        variant.INFO["{}_af".format(prefix)] = ",".join([str(self._calculate_af(global_ac[alt], global_dp)) for alt in variant.ALT])
        variant.INFO["{}_ac".format(prefix)] = ",".join([str(global_ac[alt]) for alt in variant.ALT])
        variant.INFO["{}_dp".format(prefix)] = global_dp

    def _add_insertion_stats(self, bams: List[AlignmentFile], variant: Variant, prefix: str):

        global_dp = 0
        global_ac = 0
        for i, b in enumerate(bams):
            pileups = self._get_variant_pileup(variant=variant, bam=b)
            ac, dp = get_insertion_metrics(variant=variant, pileups=pileups)
            af = self._calculate_af(ac, dp)
            if len(bams) > 1:
                variant.INFO["{}_af_{}".format(prefix, i+1)] = af
                variant.INFO["{}_ac_{}".format(prefix, i+1)] = ac
                variant.INFO["{}_dp_{}".format(prefix, i+1)] = dp
            global_ac += ac
            global_dp += dp

        global_af = self._calculate_af(global_ac, global_dp)
        variant.INFO["{}_af".format(prefix)] = global_af
        variant.INFO["{}_ac".format(prefix)] = global_ac
        variant.INFO["{}_dp".format(prefix)] = global_dp

    def _add_deletion_stats(self, bams: List[AlignmentFile], variant: Variant, prefix: str):

        global_dp = 0
        global_ac = 0
        for i, b in enumerate(bams):
            pileups = self._get_variant_pileup(variant=variant, bam=b)
            ac, dp = get_deletion_metrics(variant=variant, pileups=pileups)
            af = self._calculate_af(ac, dp)
            if len(bams) > 1:
                variant.INFO["{}_af_{}".format(prefix, i+1)] = af
                variant.INFO["{}_ac_{}".format(prefix, i+1)] = ac
                variant.INFO["{}_dp_{}".format(prefix, i+1)] = dp
            global_ac += ac
            global_dp += dp
        global_af = self._calculate_af(global_ac, global_dp)
        variant.INFO["{}_af".format(prefix)] = global_af
        variant.INFO["{}_ac".format(prefix)] = global_ac
        variant.INFO["{}_dp".format(prefix)] = global_dp

    def _calculate_af(self, ac, dp):
        return float(ac) / dp if dp > 0 else 0.0

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
