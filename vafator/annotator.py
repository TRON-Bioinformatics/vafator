from collections import Counter

import pysam
from cyvcf2 import VCF, Writer, Variant
import os
import vafator
import datetime
import json
import asyncio
import time
from vafator.power import PowerCalculator
from vafator.pileups import get_variant_pileup, get_metrics


BATCH_SIZE = 10000


def background(f):
    def wrapped(*args, **kwargs):
        return asyncio.get_event_loop().run_in_executor(None, f, *args, **kwargs)

    return wrapped


class Annotator(object):

    vafator_header = {
        "name": "vafator",
        "version": vafator.VERSION,
        "date": datetime.datetime.now().ctime(),
        "timestamp": datetime.datetime.now().timestamp(),
    }

    def __init__(self, input_vcf: str, output_vcf: str,
                 input_bams: dict,
                 purities: dict = {},
                 mapping_qual_thr=0,
                 base_call_qual_thr=29,
                 tumor_ploidies: dict = {},
                 normal_ploidy=2):

        self.mapping_quality_threshold = mapping_qual_thr
        self.base_call_quality_threshold = base_call_qual_thr
        self.purities = purities
        self.tumor_ploidies = tumor_ploidies
        self.normal_ploidy = normal_ploidy
        self.power = PowerCalculator(normal_ploidy=normal_ploidy, tumor_ploidies=tumor_ploidies, purities=purities)

        self.vcf = VCF(input_vcf)
        # sets a line in the header with the command used to annotate the file
        self.vafator_header["input_vcf"] = os.path.abspath(input_vcf)
        self.vafator_header["output_vcf"] = os.path.abspath(output_vcf)
        self.vafator_header["bams"] = ";".join(
            ["{}:{}".format(s, ",".join([os.path.abspath(b) for b in bams])) for s, bams in input_bams.items()])
        self.vafator_header["mapping_quality_threshold"] = mapping_qual_thr
        self.vafator_header["base_call_quality_threshold"] = base_call_qual_thr
        self.vafator_header["purities"] = ";".join(["{}:{}".format(s, p) for s, p in purities.items()])
        self.vafator_header["normal_ploidy"] = normal_ploidy
        self.vafator_header["tumor_ploidy"] = ";".join(["{}:{}".format(s, p) for s, p in tumor_ploidies.items()])
        self.vcf.add_to_header("##vafator_command_line={}".format(json.dumps(self.vafator_header)))
        # adds to the header all the names of the annotations
        for a in Annotator._get_headers(input_bams):
            self.vcf.add_info_to_header(a)
        self.vcf_writer = Writer(output_vcf, self.vcf)
        self.bam_readers = {s : [pysam.AlignmentFile(b, "rb") for b in bams] for s, bams in input_bams.items()}

    @staticmethod
    def _get_headers(input_bams: dict):
        headers = []

        for s, bams in input_bams.items():
            headers.append({
                'ID': "{}_af".format(s),
                'Description': "Allele frequency for the alternate alleles in the {} sample/s".format(s),
                'Type': 'Float',
                'Number': 'A'
            })
            headers.append({
                'ID': "{}_dp".format(s),
                'Description': "Total depth of coverage in the {} sample/s (independent of alleles)".format(s),
                'Type': 'Float',
                'Number': '1'
            })
            headers.append({
                'ID': "{}_ac".format(s),
                'Description': "Allele count for the alternate alleles in the {} sample/s".format(s),
                'Type': 'Integer',
                'Number': 'A'
            })
            headers.append({
                'ID': "{}_pu".format(s),
                'Description': "Probability of an undetected mutation given the observed supporting reads (AC), "
                               "the observed total coverage (DP) and the provided tumor purity in the "
                               "{} sample/s".format(s),
                'Type': 'Float',
                'Number': 'A'
            })
            headers.append({
                'ID': "{}_pw".format(s),
                'Description': "Power to detect a somatic mutation as described in Absolute "
                               "given the observed supporting reads (AC), the observed total coverage (DP) "
                               "and the provided tumor purity and ploidies in the {} sample/s".format(s),
                'Type': 'Float',
                'Number': 'A'
            })
            headers.append({
                'ID': "{}_eaf".format(s),
                'Description': "Expected VAF considering the purity and ploidy/copy number in the "
                               "{} sample/s".format(s),
                'Type': 'Float',
                'Number': '1'
            })
            headers.append({
                'ID': "{}_bq".format(s),
                'Description': "Median base call quality of the reads supporting each allele in the "
                               "{} sample/s".format(s),
                'Type': 'Float',
                'Number': 'R'
            })
            headers.append({
                'ID': "{}_mq".format(s),
                'Description': "Median mapping quality of the reads supporting each allele in the "
                               "{} sample/s".format(s),
                'Type': 'Float',
                'Number': 'R'
            })
            headers.append({
                'ID': "{}_pos".format(s),
                'Description': "Median position within the read of the reads supporting each allele in the "
                               "{} sample/s".format(s),
                'Type': 'Float',
                'Number': 'R'
            })

            if len(bams) > 1:
                for i, bam in enumerate(bams, start=1):
                    n = os.path.basename(bam).split(".")[0]
                    headers = headers + [
                        {'ID': "{}_af_{}".format(s, i),
                         'Description': "Allele frequency for the alternate alleles in the {} sample {}".format(s, n),
                         'Type': 'Float', 'Number': 'A'},
                        {'ID': "{}_dp_{}".format(s, i),
                         'Description': "Depth of coverage in the {} sample {} (independent of alleles)".format(s, n),
                         'Type': 'Float', 'Number': '1'},
                        {'ID': "{}_ac_{}".format(s, i),
                         'Description': "Allele count for the alternate alleles in the {} sample {}".format(s, n),
                         'Type': 'Integer', 'Number': 'A'},
                        {'ID': "{}_pu_{}".format(s, i),
                         'Description': "Probability of an undetected mutation given the observed supporting "
                                        "reads (AC), the observed total coverage (DP) and the provided tumor "
                                        "purity in the {} sample {}".format(s, n),
                         'Type': 'Float', 'Number': 'A'},
                        {'ID': "{}_pw_{}".format(s, i),
                         'Description': "Power to detect a somatic mutation as described in Absolute "
                                        "given the observed supporting reads (AC), the observed total coverage (DP) "
                                        "and the provided tumor purity and ploidies in the {} sample {}".format(s, n),
                         'Type': 'Float', 'Number': 'A'},
                        {'ID': "{}_bq_{}".format(s, i),
                         'Description': "Median base call quality of the reads supporting each allele in "
                                        "the {} sample {}".format(s, n),
                         'Type': 'Integer', 'Number': 'R'},
                        {'ID': "{}_mq_{}".format(s, i),
                         'Description': "Median mapping quality of the reads supporting each allele in "
                                        "the {} sample {}".format(s, n),
                         'Type': 'Integer', 'Number': 'R'},
                        {'ID': "{}_pos_{}".format(s, i),
                         'Description': "Median position within the read of the reads supporting each allele in "
                                        "the {} sample {}".format(s, n),
                         'Type': 'Integer', 'Number': 'R'},
                    ]
        return headers

    @background
    def _write_batch(self, batch):
        for v in batch:
            self.vcf_writer.write_record(v)

    def _add_stats(self, variant: Variant):
        for sample, bams in self.bam_readers.items():
            global_dp = 0
            global_ac = Counter()
            global_bq = Counter()
            global_mq = Counter()
            global_pos = Counter()
            for i, bam in enumerate(bams):
                pileups = get_variant_pileup(
                    variant=variant, bam=bam,
                    min_base_quality=self.base_call_quality_threshold,
                    min_mapping_quality=self.mapping_quality_threshold)
                coverage_metrics = get_metrics(variant=variant, pileups=pileups)
                if coverage_metrics is not None:
                    if len(bams) > 1:
                        variant.INFO["{}_af_{}".format(sample, i + 1)] = ",".join(
                            [str(self._calculate_af(coverage_metrics.ac[alt], coverage_metrics.dp)) for alt in variant.ALT])
                        variant.INFO["{}_ac_{}".format(sample, i + 1)] = ",".join([str(coverage_metrics.ac[alt]) for alt in variant.ALT])
                        variant.INFO["{}_dp_{}".format(sample, i + 1)] = coverage_metrics.dp
                        variant.INFO["{}_pu_{}".format(sample, i + 1)] = ",".join(
                            [str(self.power.calculate_power(
                                ac=coverage_metrics.ac[alt], dp=coverage_metrics.dp, sample=sample, variant=variant
                            )) for alt in variant.ALT])
                        variant.INFO["{}_pw_{}".format(sample, i + 1)] = ",".join(
                            [str(self.power.calculate_absolute_power(
                                ac=coverage_metrics.ac[alt], dp=coverage_metrics.dp, sample=sample, variant=variant
                            )) for alt in variant.ALT])
                        variant.INFO["{}_bq_{}".format(sample, i + 1)] = ",".join(
                            [str(coverage_metrics.bqs[variant.REF])] +
                            [str(coverage_metrics.bqs[alt]) for alt in variant.ALT])
                        variant.INFO["{}_mq_{}".format(sample, i + 1)] = ",".join(
                            [str(coverage_metrics.mqs[variant.REF])] +
                            [str(coverage_metrics.mqs[alt]) for alt in variant.ALT])
                        variant.INFO["{}_pos_{}".format(sample, i + 1)] = ",".join(
                            [str(coverage_metrics.positions[variant.REF])] +
                            [str(coverage_metrics.positions[alt]) for alt in variant.ALT])
                    global_ac.update(coverage_metrics.ac)
                    global_bq.update(coverage_metrics.bqs)
                    global_mq.update(coverage_metrics.mqs)
                    global_pos.update(coverage_metrics.positions)
                    global_dp += coverage_metrics.dp

            variant.INFO["{}_af".format(sample)] = ",".join([str(self._calculate_af(global_ac[alt], global_dp)) for alt in variant.ALT])
            variant.INFO["{}_ac".format(sample)] = ",".join([str(global_ac[alt]) for alt in variant.ALT])
            variant.INFO["{}_dp".format(sample)] = global_dp
            variant.INFO["{}_eaf".format(sample)] = str(self.power.calculate_expected_vaf(
                sample=pysam, variant=variant))
            variant.INFO["{}_pu".format(sample)] = ",".join(
                [str(self.power.calculate_power(ac=global_ac[alt], dp=global_dp, sample=sample, variant=variant))
                 for alt in variant.ALT])
            variant.INFO["{}_pw".format(sample)] = ",".join(
                [str(self.power.calculate_absolute_power(
                    ac=global_ac[alt], dp=global_dp, sample=sample, variant=variant))
                 for alt in variant.ALT])
            variant.INFO["{}_bq".format(sample)] = ",".join(
                [str(global_bq[variant.REF])] + [str(global_bq[alt]) for alt in variant.ALT])
            variant.INFO["{}_mq".format(sample)] = ",".join(
                [str(global_mq[variant.REF])] + [str(global_mq[alt]) for alt in variant.ALT])
            variant.INFO["{}_pos".format(sample)] = ",".join(
                [str(global_pos[variant.REF])] + [str(global_pos[alt]) for alt in variant.ALT])

    def _calculate_af(self, ac, dp):
        return round(float(ac) / dp, 5) if dp > 0 else 0.0

    def run(self):
        batch = []
        variant: Variant
        for variant in self.vcf:
            # gets the counts of all bases across all BAMs
            self._add_stats(variant)

            batch.append(variant)
            if len(batch) >= BATCH_SIZE:
                self._write_batch(batch)
                batch = []
        if len(batch) > 0:
            self._write_batch(batch)

        time.sleep(2)

        self.vcf_writer.close()
        self.vcf.close()
        for _, bams in self.bam_readers.items():
            for bam in bams:
                bam.close()
