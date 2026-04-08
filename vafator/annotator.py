from collections import Counter
import os
from concurrent.futures import ProcessPoolExecutor

import pysam
from cyvcf2 import VCF, Writer, Variant
import vafator
import datetime
import json

from vafator.constants import (
    AMBIGUOUS_BASES,
    BATCH_SIZE,
    _HEADER_TEMPLATES,
    _REPLICATE_HEADER_TEMPLATES,
)
from vafator.ploidies import DEFAULT_PLOIDY
from vafator.rank_sum_test import get_rank_sum_tests
from vafator.power import PowerCalculator, DEFAULT_ERROR_RATE, DEFAULT_FPR
from vafator.pileups import collect_metrics_for_chrom, stream_variants_by_chrom
from vafator.pileup_utils import VariantRecord, EMPTY_METRICS


def _collect_metrics_worker(
    chrom: str,
    variant_tuples: list,
    bam_paths: dict,
    min_base_quality: int,
    min_mapping_quality: int,
    include_ambiguous_bases: bool,
) -> dict:
    """
    Top-level worker function for ProcessPoolExecutor — must be module-level to be picklable.
    Opens its own BAM readers (AlignmentFile objects cannot be shared across processes).
    Receives variant data as plain tuples (cyvcf2.Variant objects are not picklable).

    Args:
        chrom: chromosome name
        variant_tuples: list of (POS, REF, ALT) tuples for variants on this chromosome
        bam_paths: {sample: [bam_path, ...]} — picklable BAM file paths
        min_base_quality: minimum base call quality threshold
        min_mapping_quality: minimum mapping quality threshold
        include_ambiguous_bases: whether to include ambiguous bases in depth calculation

    Returns:
        {(sample, bam_idx): {(pos, REF, ALT): CoverageMetrics}}
    """
    all_metrics = {}
    variants = [
        VariantRecord(CHROM=chrom, POS=pos, REF=ref, ALT=[alt])
        for pos, ref, alt in variant_tuples
    ]
    for sample, bam_files in bam_paths.items():
        for i, bam_path in enumerate(bam_files):
            bam = pysam.AlignmentFile(bam_path, "rb")
            all_metrics[(sample, i)] = collect_metrics_for_chrom(
                chrom=chrom,
                variants=variants,
                bam=bam,
                min_base_quality=min_base_quality,
                min_mapping_quality=min_mapping_quality,
                include_ambiguous_bases=include_ambiguous_bases,
            )
            bam.close()
    return all_metrics


class Annotator(object):

    def __init__(
        self,
        input_vcf: str,
        output_vcf: str,
        input_bams: dict,
        purities: dict = {},
        mapping_qual_thr: int = 0,
        base_call_qual_thr: int = 29,
        tumor_ploidies: dict = {},
        normal_ploidy: int = 2,
        fpr: float = DEFAULT_FPR,
        error_rate: float = DEFAULT_ERROR_RATE,
        include_ambiguous_bases: bool = True,
        num_processes: int = 1,
    ):
        """
        Args:
            input_vcf: path to the input VCF file to annotate
            output_vcf: path for the annotated output VCF
            input_bams: {sample_name: [bam_path, ...]} — one or more BAMs per sample
            purities: {sample_name: purity} — tumor purity per sample (default: 1.0)
            mapping_qual_thr: minimum mapping quality; reads below this are excluded
            base_call_qual_thr: minimum base call quality; bases below this are excluded
            tumor_ploidies: {sample_name: PloidyManager} — tumor ploidy per sample (default: 2)
            normal_ploidy: normal ploidy for power calculation (default: 2)
            fpr: false positive rate for power calculation
            error_rate: sequencing error rate for power calculation
            include_ambiguous_bases: if True, ambiguous bases (N and IUPAC codes) are counted in DP
            num_processes: number of parallel processes for chromosome-level annotation (default: 1)
        """

        self.mapping_quality_threshold = mapping_qual_thr
        self.base_call_quality_threshold = base_call_qual_thr
        self.purities = purities
        self.tumor_ploidies = tumor_ploidies
        self.normal_ploidy = normal_ploidy
        self.include_ambiguous_bases = include_ambiguous_bases
        self.num_processes = num_processes
        self.power = PowerCalculator(
            normal_ploidy=normal_ploidy,
            tumor_ploidies=tumor_ploidies,
            purities=purities,
            error_rate=error_rate,
            fpr=fpr,
        )

        self.vcf = VCF(input_vcf)

        self.header = {
            "name": "vafator",
            "version": vafator.VERSION,
            "date": datetime.datetime.now().ctime(),
            "timestamp": datetime.datetime.now().timestamp(),
            "input_vcf": os.path.abspath(input_vcf),
            "output_vcf": os.path.abspath(output_vcf),
            "bams": ";".join(
                [
                    "{}:{}".format(s, ",".join([os.path.abspath(b) for b in bams]))
                    for s, bams in input_bams.items()
                ]
            ),
            "mapping_quality_threshold": mapping_qual_thr,
            "base_call_quality_threshold": base_call_qual_thr,
            "purities": ";".join(["{}:{}".format(s, p) for s, p in purities.items()]),
            "normal_ploidy": normal_ploidy,
            "tumor_ploidy": (
                ";".join(
                    [
                        "{}:{}".format(s, p.report_value)
                        for s, p in tumor_ploidies.items()
                    ]
                )
                if tumor_ploidies
                else DEFAULT_PLOIDY
            ),
            "include_ambiguous_bases": include_ambiguous_bases,
        }
        self.vcf.add_to_header(
            "##vafator_command_line={}".format(json.dumps(self.header))
        )

        for a in Annotator._get_headers(input_bams):
            self.vcf.add_info_to_header(a)
        self.vcf_writer = Writer(output_vcf, self.vcf)

        self.bam_paths = input_bams
        self.bam_readers = {
            s: [pysam.AlignmentFile(b, "rb") for b in bams]
            for s, bams in input_bams.items()
        }

    def run(self) -> None:
        """Run the annotation pipeline over all variants in the input VCF,
        writing annotated records to the output VCF."""
        batch = []
        if self.num_processes > 1:
            self._run_parallel(batch)
        else:
            self._run_serial(batch)
        if batch:
            self._write_batch(batch)
        self.vcf_writer.close()
        self.vcf.close()
        for _, bams in self.bam_readers.items():
            for bam in bams:
                bam.close()

    def _run_serial(self, batch: list) -> None:
        """Annotate variants chromosome by chromosome in the main process."""
        for chrom, chrom_variants in stream_variants_by_chrom(self.vcf):
            all_metrics = self._collect_chrom_metrics(chrom, chrom_variants)
            self._annotate_and_batch(chrom_variants, all_metrics, batch)

    def _run_parallel(self, batch: list) -> None:
        """Annotate variants using a process pool — one worker per chromosome.
        Results are collected and written in original VCF chromosome order."""
        chrom_variants_map = {}
        futures = {}
        with ProcessPoolExecutor(max_workers=self.num_processes) as executor:
            for chrom, chrom_variants in stream_variants_by_chrom(self.vcf):
                chrom_variants_map[chrom] = chrom_variants
                variant_tuples = [(v.POS, v.REF, v.ALT[0]) for v in chrom_variants]
                futures[
                    executor.submit(
                        _collect_metrics_worker,
                        chrom=chrom,
                        variant_tuples=variant_tuples,
                        bam_paths=self.bam_paths,
                        min_base_quality=self.base_call_quality_threshold,
                        min_mapping_quality=self.mapping_quality_threshold,
                        include_ambiguous_bases=self.include_ambiguous_bases,
                    )
                ] = chrom
            chrom_results = {futures[f]: f.result() for f in futures}
        for chrom, chrom_variants in chrom_variants_map.items():
            self._annotate_and_batch(chrom_variants, chrom_results[chrom], batch)

    def _collect_chrom_metrics(self, chrom: str, chrom_variants: list) -> dict:
        """Collect metrics for all BAMs for one chromosome in the main process.

        Args:
            chrom: chromosome name
            chrom_variants: list of cyvcf2 Variant objects on this chromosome

        Returns:
            {(sample, bam_idx): {(pos, REF, ALT): CoverageMetrics}}
        """
        all_metrics = {}
        for sample, bams in self.bam_readers.items():
            for i, bam in enumerate(bams):
                all_metrics[(sample, i)] = collect_metrics_for_chrom(
                    chrom=chrom,
                    variants=chrom_variants,
                    bam=bam,
                    min_base_quality=self.base_call_quality_threshold,
                    min_mapping_quality=self.mapping_quality_threshold,
                    include_ambiguous_bases=self.include_ambiguous_bases,
                )
        return all_metrics

    def _annotate_and_batch(
        self, chrom_variants: list, all_metrics: dict, batch: list
    ) -> None:
        """Annotate variants using pre-computed metrics and append to write batch.
        Flushes the batch to disk when it reaches BATCH_SIZE.

        Args:
            chrom_variants: list of cyvcf2 Variant objects to annotate
            all_metrics: {(sample, bam_idx): {(pos, REF, ALT): CoverageMetrics}}
            batch: accumulator list for annotated variants pending write
        """
        for variant in chrom_variants:
            metrics_by_bam = {
                (sample, i): all_metrics[(sample, i)].get(
                    (variant.POS, variant.REF, variant.ALT[0]), EMPTY_METRICS
                )
                for sample, bams in self.bam_readers.items()
                for i in range(len(bams))
            }
            self._add_stats(variant, metrics_by_bam)
            batch.append(variant)
            if len(batch) >= BATCH_SIZE:
                self._write_batch(batch)
                batch.clear()

    def _add_stats(self, variant: Variant, metrics_by_bam: dict) -> None:
        """Annotate a single variant using pre-computed metrics.

        Args:
            variant: the cyvcf2 Variant to annotate in place
            metrics_by_bam: {(sample, bam_idx): CoverageMetrics}
        """
        for sample, bams in self.bam_readers.items():
            global_dp = 0
            global_ac = Counter()
            global_bq = Counter()
            global_mq = Counter()
            global_pos = Counter()
            global_all_mqs = {}
            global_all_bqs = {}
            global_all_positions = {}

            for i in range(len(bams)):
                metrics = metrics_by_bam.get((sample, i), EMPTY_METRICS)
                if metrics is not None:
                    if len(bams) > 1:
                        self._annotate_replicate(variant, sample, i, metrics)
                    global_ac.update(metrics.ac)
                    global_bq.update(metrics.bqs)
                    global_mq.update(metrics.mqs)
                    global_pos.update(metrics.positions)
                    global_all_mqs.update(metrics.all_mqs)
                    global_all_bqs.update(metrics.all_bqs)
                    global_all_positions.update(metrics.all_positions)
                    global_dp += metrics.dp

            self._annotate_sample(
                variant,
                sample,
                global_ac,
                global_dp,
                global_bq,
                global_mq,
                global_pos,
                global_all_mqs,
                global_all_bqs,
                global_all_positions,
            )

    def _annotate_replicate(self, v: Variant, s: str, i: int, m) -> None:
        """Write per-replicate annotations — only called when multiple BAMs are provided for a sample.

        Args:
            v: the cyvcf2 Variant being annotated
            s: sample name (e.g. 'tumor', 'normal')
            i: 0-based index of the BAM within the sample's BAM list
            m: pre-computed CoverageMetrics for this variant in this BAM
        """
        n = i + 1
        v.INFO["{}_af_{}".format(s, n)] = ",".join(
            [str(self._calculate_af(m.ac[alt], m.dp)) for alt in v.ALT]
        )
        v.INFO["{}_ac_{}".format(s, n)] = ",".join([str(m.ac[alt]) for alt in v.ALT])
        v.INFO["{}_n_{}".format(s, n)] = str(
            sum(m.ac.get(b, 0) for b in AMBIGUOUS_BASES)
        )
        v.INFO["{}_dp_{}".format(s, n)] = m.dp
        v.INFO["{}_pu_{}".format(s, n)] = ",".join(
            [
                str(
                    self.power.calculate_power(
                        ac=m.ac[alt], dp=m.dp, sample=s, variant=v
                    )
                )
                for alt in v.ALT
            ]
        )

        power, k = self.power.calculate_absolute_power(dp=m.dp, sample=s, variant=v)
        v.INFO["{}_pw_{}".format(s, n)] = str(power)
        v.INFO["{}_k_{}".format(s, n)] = str(k)
        v.INFO["{}_bq_{}".format(s, n)] = ",".join(
            [str(m.bqs[v.REF])] + [str(m.bqs[alt]) for alt in v.ALT]
        )
        v.INFO["{}_mq_{}".format(s, n)] = ",".join(
            [str(m.mqs[v.REF])] + [str(m.mqs[alt]) for alt in v.ALT]
        )
        v.INFO["{}_pos_{}".format(s, n)] = ",".join(
            [str(m.positions[v.REF])] + [str(m.positions[alt]) for alt in v.ALT]
        )
        for key, tag in [
            (m.all_mqs, "rsmq"),
            (m.all_bqs, "rsbq"),
            (m.all_positions, "rspos"),
        ]:
            pvalues, stats = get_rank_sum_tests(key, v)
            if stats:
                v.INFO["{}_{}_{}".format(s, tag, n)] = ",".join(stats)
                v.INFO["{}_{}_pv_{}".format(s, tag, n)] = ",".join(pvalues)

    def _annotate_sample(
        self,
        v: Variant,
        s: str,
        gac: Counter,
        gdp: int,
        gbq: Counter,
        gmq: Counter,
        gpos: Counter,
        gallmq: dict,
        gallbq: dict,
        gallpos: dict,
    ) -> None:
        """Write aggregate annotations for a sample, combining metrics across all replicates.

        Args:
            v: the cyvcf2 Variant being annotated
            s: sample name (e.g. 'tumor', 'normal')
            gac: allele counts summed across all BAMs
            gdp: total depth summed across all BAMs
            gbq: median BQ per allele, summed across all BAMs
            gmq: median MQ per allele, summed across all BAMs
            gpos: median read position per allele, summed across all BAMs
            gallmq: MQ distributions per allele across all BAMs
            gallbq: BQ distributions per allele across all BAMs
            gallpos: read position distributions per allele across all BAMs
        """
        v.INFO["{}_af".format(s)] = ",".join(
            [str(self._calculate_af(gac[alt], gdp)) for alt in v.ALT]
        )
        v.INFO["{}_ac".format(s)] = ",".join([str(gac[alt]) for alt in v.ALT])
        v.INFO["{}_n".format(s)] = str(sum(gac.get(b, 0) for b in AMBIGUOUS_BASES))
        v.INFO["{}_dp".format(s)] = gdp
        v.INFO["{}_eaf".format(s)] = str(
            self.power.calculate_expected_vaf(sample=s, variant=v)
        )
        v.INFO["{}_pu".format(s)] = ",".join(
            [
                str(
                    self.power.calculate_power(ac=gac[alt], dp=gdp, sample=s, variant=v)
                )
                for alt in v.ALT
            ]
        )

        power, k = self.power.calculate_absolute_power(dp=gdp, sample=s, variant=v)
        v.INFO["{}_pw".format(s)] = str(power)
        v.INFO["{}_k".format(s)] = str(k)
        v.INFO["{}_bq".format(s)] = ",".join(
            [str(gbq[v.REF])] + [str(gbq[alt]) for alt in v.ALT]
        )
        v.INFO["{}_mq".format(s)] = ",".join(
            [str(gmq[v.REF])] + [str(gmq[alt]) for alt in v.ALT]
        )
        v.INFO["{}_pos".format(s)] = ",".join(
            [str(gpos[v.REF])] + [str(gpos[alt]) for alt in v.ALT]
        )
        for distributions, tag in [
            (gallmq, "rsmq"),
            (gallbq, "rsbq"),
            (gallpos, "rspos"),
        ]:
            pvalues, stats = get_rank_sum_tests(distributions, v)
            if stats:
                v.INFO["{}_{}".format(s, tag)] = ",".join(stats)
                v.INFO["{}_{}_pv".format(s, tag)] = ",".join(pvalues)

    def _calculate_af(self, ac: int, dp: int) -> float:
        """Return allele frequency, or 0.0 if depth is zero.

        Args:
            ac: allele count for this alternate allele
            dp: total depth of coverage

        Returns:
            allele frequency rounded to 5 decimal places, or 0.0 if dp is zero
        """
        return round(float(ac) / dp, 5) if dp > 0 else 0.0

    def _write_batch(self, batch: list) -> None:
        """Write a batch of annotated variants to the output VCF.

        Args:
            batch: list of cyvcf2 Variant objects to write
        """
        for v in batch:
            self.vcf_writer.write_record(v)

    @staticmethod
    def _get_headers(input_bams: dict) -> list:
        """Build the list of INFO header entries for all samples and replicates.

        Args:
            input_bams: {sample_name: [bam_path, ...]}

        Returns:
            list of header dicts suitable for cyvcf2's add_info_to_header
        """
        headers = []
        for s, bams in input_bams.items():
            for suffix, description, typ, number in _HEADER_TEMPLATES:
                headers.append(
                    Annotator._make_header(suffix, description, typ, number, sample=s)
                )
            if len(bams) > 1:
                for i in range(1, len(bams) + 1):
                    # n = os.path.basename(bam).split(".")[0]
                    for suffix, description, typ, number in _REPLICATE_HEADER_TEMPLATES:
                        headers.append(
                            Annotator._make_header(
                                suffix, description, typ, number, sample=s, index=i
                            )
                        )
        return headers

    @staticmethod
    def _make_header(
        suffix: str,
        description: str,
        typ: str,
        number: str,
        sample: str,
        index: int = None,
    ) -> dict:
        """Build a single INFO header dict for cyvcf2's add_info_to_header.

        Args:
            suffix: annotation suffix (e.g. 'af', 'dp')
            description: description template with {sample} placeholder
            typ: VCF type string ('Float', 'Integer', 'String')
            number: VCF number string ('A', 'R', '1', etc.)
            sample: sample name to substitute into ID and description
            index: optional replicate index appended to the ID

        Returns:
            dict with keys ID, Description, Type, Number suitable for cyvcf2's add_info_to_header
        """
        header_id = "{}_{}".format(sample, suffix)
        if index is not None:
            header_id = "{}_{}".format(header_id, index)
        return {
            "ID": header_id,
            "Description": description.format(sample=sample),
            "Type": typ,
            "Number": number,
        }
