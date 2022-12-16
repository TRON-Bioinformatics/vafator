from collections import Counter

import numpy as np
import pysam
from cyvcf2 import VCF, Writer, Variant
import os
import vafator
import datetime
import json
import asyncio
import time

from vafator.ploidies import DEFAULT_PLOIDY
from vafator.rank_sum_test import calculate_rank_sum_test, get_rank_sum_tests
from vafator.power import PowerCalculator, DEFAULT_ERROR_RATE, DEFAULT_FPR
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
                 normal_ploidy=2,
                 fpr=DEFAULT_FPR,
                 error_rate=DEFAULT_ERROR_RATE,
                 include_ambiguous_bases=False):

        self.mapping_quality_threshold = mapping_qual_thr
        self.base_call_quality_threshold = base_call_qual_thr
        self.purities = purities
        self.tumor_ploidies = tumor_ploidies
        self.normal_ploidy = normal_ploidy
        self.include_ambiguous_bases = include_ambiguous_bases
        self.power = PowerCalculator(
            normal_ploidy=normal_ploidy, tumor_ploidies=tumor_ploidies, purities=purities,
            error_rate=error_rate, fpr=fpr)

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
        self.vafator_header["tumor_ploidy"] = ";".join(["{}:{}".format(s, p.report_value)
                                                        for s, p in tumor_ploidies.items()]) \
            if tumor_ploidies else DEFAULT_PLOIDY
        self.vafator_header["include_ambiguous_bases"] = self.include_ambiguous_bases
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
                'ID': "{}_n".format(s),
                'Description': "Allele count for ambiguous bases (any IUPAC ambiguity code is counted) "
                               "in the {} sample/s".format(s),
                'Type': 'Integer',
                'Number': '1'
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
                               "given the observed total coverage (DP) "
                               "and the provided tumor purity and ploidies in the {} sample/s".format(s),
                'Type': 'Float',
                'Number': '1'
            })
            headers.append({
                'ID': "{}_k".format(s),
                'Description': "Minimum number of supporting reads, k, such that the probability of observing "
                               "k or more non-reference reads due to sequencing error is less than the defined FPR "
                               "in the {} sample/s".format(s),
                'Type': 'Float',
                'Number': '1'
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
            headers.append({
                'ID': "{}_rsmq".format(s),
                'Description': "Rank sum test comparing the MQ distributions supporting the reference and the "
                               "alternate in the {} sample/s. Identical distributions will have a value of 0, larger "
                               "values away from 0 indicate different distributions.".format(s),
                'Type': 'Float',
                'Number': 'A'
            })
            headers.append({
                'ID': "{}_rsmq_pv".format(s),
                'Description': "Rank sum test comparing the mapping quality distributions between alternate "
                               "and reference p-value in the {} sample/s. , The null hypothesis is that there is no "
                               "difference between the distributions".format(s),
                'Type': 'Float',
                'Number': 'A'
            })
            headers.append({
                'ID': "{}_rsbq".format(s),
                'Description': "Rank sum test comparing the base call qualities distributions supporting the reference "
                               "and the alternate in the {} sample/s. Identical distributions will have a value of 0, "
                               "larger values away from 0 indicate different distributions.".format(s),
                'Type': 'Float',
                'Number': 'A'
            })
            headers.append({
                'ID': "{}_rsbq_pv".format(s),
                'Description': "Rank sum test comparing the base call qualities distributions between alternate "
                               "and reference p-value in the {} sample/s. , The null hypothesis is that there is no "
                               "difference between the distributions".format(s),
                'Type': 'Float',
                'Number': 'A'
            })
            headers.append({
                'ID': "{}_rspos".format(s),
                'Description': "Rank sum test comparing the relative position distributions supporting the reference "
                               "and the alternate in the {} sample/s. Identical distributions will have a value of 0, "
                               "larger values away from 0 indicate different distributions.".format(s),
                'Type': 'Float',
                'Number': 'A'
            })
            headers.append({
                'ID': "{}_rspos_pv".format(s),
                'Description': "Rank sum test comparing the relative position distributions between alternate "
                               "and reference p-value in the {} sample/s. , The null hypothesis is that there is no "
                               "difference between the distributions".format(s),
                'Type': 'Float',
                'Number': 'A'
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
                        {'ID': "{}_n_{}".format(s, i),
                         'Description': "Allele count for ambiguous bases (any IUPAC ambiguity code is counted) "
                                        "in the {} sample {}".format(s, n),
                         'Type': 'Integer', 'Number': '1'},
                        {'ID': "{}_pu_{}".format(s, i),
                         'Description': "Probability of an undetected mutation given the observed supporting "
                                        "reads (AC), the observed total coverage (DP) and the provided tumor "
                                        "purity in the {} sample {}".format(s, n),
                         'Type': 'Float', 'Number': 'A'},
                        {'ID': "{}_pw_{}".format(s, i),
                         'Description': "Power to detect a somatic mutation as described in Absolute "
                                        "given the observed total coverage (DP) "
                                        "and the provided tumor purity and ploidies in the {} sample {}".format(s, n),
                         'Type': 'Float', 'Number': '1'},
                        {'ID': "{}_k_{}".format(s, i),
                         'Description': "Minimum number of supporting reads, k, such that the probability of observing "
                                        "k or more non-reference reads due to sequencing error is less than the "
                                        "defined FPR in the {} sample {}".format(s, n),
                         'Type': 'Float',
                         'Number': '1'},
                        {'ID': "{}_bq_{}".format(s, i),
                         'Description': "Median base call quality of the reads supporting each allele in "
                                        "the {} sample {}".format(s, n),
                         'Type': 'Float', 'Number': 'R'},
                        {'ID': "{}_rsbq_{}".format(s, i),
                         'Description': "Rank sum test comparing the base call qualities distributions supporting the "
                                        "reference and the alternate in the {} sample {}".format(s, n),
                         'Type': 'Float', 'Number': 'A'},
                        {'ID': "{}_rsbq_pv_{}".format(s, i),
                         'Description': "Significance for the rank sum test comparing the base call qualities "
                                        "distributions supporting the reference and the alternate "
                                        "in the {} sample {}".format(s, n),
                         'Type': 'Float', 'Number': 'A'},
                        {'ID': "{}_mq_{}".format(s, i),
                         'Description': "Median mapping quality of the reads supporting each allele in "
                                        "the {} sample {}".format(s, n),
                         'Type': 'Float', 'Number': 'R'},
                        {'ID': "{}_rsmq_{}".format(s, i),
                         'Description': "Rank sum test comparing the mapping qualities distributions supporting the "
                                        "reference and the alternate in the {} sample {}".format(s, n),
                         'Type': 'Float', 'Number': 'A'},
                        {'ID': "{}_rsmq_pv_{}".format(s, i),
                         'Description': "Significance for the rank sum test comparing the mapping qualities "
                                        "distributions supporting the reference and the alternate "
                                        "in the {} sample {}".format(s, n),
                         'Type': 'Float', 'Number': 'A'},
                        {'ID': "{}_pos_{}".format(s, i),
                         'Description': "Median position within the read of the reads supporting each allele in "
                                        "the {} sample {}".format(s, n),
                         'Type': 'Float', 'Number': 'R'},
                        {'ID': "{}_rspos_{}".format(s, i),
                         'Description': "Rank sum test comparing the position distributions supporting the "
                                        "reference and the alternate in the {} sample {}".format(s, n),
                         'Type': 'Float', 'Number': 'A'},
                        {'ID': "{}_rspos_pv_{}".format(s, i),
                         'Description': "Significance for the rank sum test comparing the position "
                                        "distributions supporting the reference and the alternate "
                                        "in the {} sample {}".format(s, n),
                         'Type': 'Float', 'Number': 'A'},
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
            global_all_mqs = {}
            global_all_bqs = {}
            global_all_positions = {}
            for i, bam in enumerate(bams):
                pileups = get_variant_pileup(
                    variant=variant, bam=bam,
                    min_base_quality=self.base_call_quality_threshold,
                    min_mapping_quality=self.mapping_quality_threshold)
                coverage_metrics = get_metrics(variant=variant, pileups=pileups,
                                               include_ambiguous_bases=self.include_ambiguous_bases)
                if coverage_metrics is not None:
                    if len(bams) > 1:
                        variant.INFO["{}_af_{}".format(sample, i + 1)] = \
                            ",".join([str(self._calculate_af(coverage_metrics.ac[alt], coverage_metrics.dp))
                                      for alt in variant.ALT])
                        variant.INFO["{}_ac_{}".format(sample, i + 1)] = \
                            ",".join([str(coverage_metrics.ac[alt]) for alt in variant.ALT])
                        variant.INFO["{}_n_{}".format(sample, i + 1)] = \
                            str(sum([coverage_metrics.ac.get(n, 0) for n in vafator.AMBIGUOUS_BASES]))
                        variant.INFO["{}_dp_{}".format(sample, i + 1)] = coverage_metrics.dp
                        variant.INFO["{}_pu_{}".format(sample, i + 1)] = ",".join(
                            [str(self.power.calculate_power(
                                ac=coverage_metrics.ac[alt], dp=coverage_metrics.dp, sample=sample, variant=variant
                            )) for alt in variant.ALT])
                        power, k = self.power.calculate_absolute_power(
                            dp=coverage_metrics.dp, sample=sample, variant=variant)
                        variant.INFO["{}_pw_{}".format(sample, i + 1)] = str(power)
                        variant.INFO["{}_k_{}".format(sample, i + 1)] = str(k)
                        variant.INFO["{}_bq_{}".format(sample, i + 1)] = ",".join(
                            [str(coverage_metrics.bqs[variant.REF])] +
                            [str(coverage_metrics.bqs[alt]) for alt in variant.ALT])
                        variant.INFO["{}_mq_{}".format(sample, i + 1)] = ",".join(
                            [str(coverage_metrics.mqs[variant.REF])] +
                            [str(coverage_metrics.mqs[alt]) for alt in variant.ALT])
                        variant.INFO["{}_pos_{}".format(sample, i + 1)] = ",".join(
                            [str(coverage_metrics.positions[variant.REF])] +
                            [str(coverage_metrics.positions[alt]) for alt in variant.ALT])

                        pvalues, stats = get_rank_sum_tests(coverage_metrics.all_mqs, variant)
                        if stats:
                            variant.INFO["{}_rsmq_{}".format(sample, i + 1)] = ",".join(stats)
                            variant.INFO["{}_rsmq_pv_{}".format(sample, i + 1)] = ",".join(pvalues)

                        pvalues, stats = get_rank_sum_tests(coverage_metrics.all_bqs, variant)
                        if stats:
                            variant.INFO["{}_rsbq_{}".format(sample, i + 1)] = ",".join(stats)
                            variant.INFO["{}_rsbq_pv_{}".format(sample, i + 1)] = ",".join(pvalues)

                        pvalues, stats = get_rank_sum_tests(coverage_metrics.all_positions, variant)
                        if stats:
                            variant.INFO["{}_rspos_{}".format(sample, i + 1)] = ",".join(stats)
                            variant.INFO["{}_rspos_pv_{}".format(sample, i + 1)] = ",".join(pvalues)

                    global_ac.update(coverage_metrics.ac)
                    global_bq.update(coverage_metrics.bqs)
                    global_mq.update(coverage_metrics.mqs)
                    global_pos.update(coverage_metrics.positions)
                    global_all_mqs.update(coverage_metrics.all_mqs)
                    global_all_bqs.update(coverage_metrics.all_bqs)
                    global_all_positions.update(coverage_metrics.all_positions)
                    global_dp += coverage_metrics.dp

            variant.INFO["{}_af".format(sample)] = ",".join([str(self._calculate_af(global_ac[alt], global_dp)) for alt in variant.ALT])
            variant.INFO["{}_ac".format(sample)] = ",".join([str(global_ac[alt]) for alt in variant.ALT])
            variant.INFO["{}_n".format(sample)] = str(sum([global_ac.get(n, 0) for n in vafator.AMBIGUOUS_BASES]))
            variant.INFO["{}_dp".format(sample)] = global_dp
            variant.INFO["{}_eaf".format(sample)] = str(self.power.calculate_expected_vaf(
                sample=sample, variant=variant))
            variant.INFO["{}_pu".format(sample)] = ",".join(
                [str(self.power.calculate_power(ac=global_ac[alt], dp=global_dp, sample=sample, variant=variant))
                 for alt in variant.ALT])
            power, k = self.power.calculate_absolute_power(
                dp=global_dp, sample=sample, variant=variant)
            variant.INFO["{}_pw".format(sample)] = str(power)
            variant.INFO["{}_k".format(sample)] = str(k)
            variant.INFO["{}_bq".format(sample)] = ",".join(
                [str(global_bq[variant.REF])] + [str(global_bq[alt]) for alt in variant.ALT])
            variant.INFO["{}_mq".format(sample)] = ",".join(
                [str(global_mq[variant.REF])] + [str(global_mq[alt]) for alt in variant.ALT])
            variant.INFO["{}_pos".format(sample)] = ",".join(
                [str(global_pos[variant.REF])] + [str(global_pos[alt]) for alt in variant.ALT])

            # for these rank sum tests it is required at least one value for the alternate and one value for the
            # reference otherwise it cannot be calculated
            pvalues, stats = get_rank_sum_tests(global_all_mqs, variant)
            if stats:
                variant.INFO["{}_rsmq".format(sample)] = ",".join(stats)
                variant.INFO["{}_rsmq_pv".format(sample)] = ",".join(pvalues)

            pvalues, stats = get_rank_sum_tests(global_all_bqs, variant)
            if stats:
                variant.INFO["{}_rsbq".format(sample)] = ",".join(stats)
                variant.INFO["{}_rsbq_pv".format(sample)] = ",".join(pvalues)

            pvalues, stats = get_rank_sum_tests(global_all_positions, variant)
            if stats:
                variant.INFO["{}_rspos".format(sample)] = ",".join(stats)
                variant.INFO["{}_rspos_pv".format(sample)] = ",".join(pvalues)

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
