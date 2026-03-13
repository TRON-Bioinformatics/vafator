from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Union, List, Dict, Iterator, Tuple
from cyvcf2 import Variant
from pysam.libcalignmentfile import IteratorColumnRegion, AlignmentFile

from vafator import AMBIGUOUS_BASES
from vafator.tests.utils import VafatorVariant
import numpy as np


def is_snp(variant: Variant):
    return len(variant.REF) == 1 and len(variant.ALT[0]) == 1


def is_insertion(variant: Variant):
    return len(variant.REF) == 1 and len(variant.ALT[0]) > 1


def is_deletion(variant: Variant):
    return len(variant.ALT[0]) == 1 and len(variant.REF) > 1


def get_variant_pileup(
        variant: Union[Variant, VafatorVariant], bam: AlignmentFile,
        min_base_quality, min_mapping_quality) -> IteratorColumnRegion:
    """Single-variant pileup, kept for backwards compatibility and tests."""
    position = variant.POS
    return bam.pileup(contig=variant.CHROM, start=position - 1, stop=position,
                      truncate=True,
                      max_depth=1000000,
                      min_base_quality=min_base_quality,
                      min_mapping_quality=min_mapping_quality,
                      stepper='samtools',
                      )


def get_region_pileup(chrom: str, start: int, end: int, bam: AlignmentFile,
                      min_base_quality, min_mapping_quality):
    """
    Opens a single pileup iterator spanning a whole region (e.g. one chromosome).
    start is 0-based inclusive, end is 1-based exclusive (last variant POS).
    """
    return bam.pileup(contig=chrom, start=start, stop=end,
                      truncate=True,
                      max_depth=1000000,
                      min_base_quality=min_base_quality,
                      min_mapping_quality=min_mapping_quality,
                      stepper='samtools',
                      )


def stream_variants_by_chrom(vcf) -> Iterator[Tuple[str, List[Variant]]]:
    """
    Yields (chrom, [variants]) one chromosome at a time.
    Only one chromosome's variants are held in memory at once.
    """
    current_chrom = None
    current_variants = []
    for variant in vcf:
        if variant.CHROM != current_chrom:
            if current_variants:
                yield current_chrom, current_variants
            current_chrom = variant.CHROM
            current_variants = [variant]
        else:
            current_variants.append(variant)
    if current_variants:
        yield current_chrom, current_variants


def collect_metrics_for_chrom(
        chrom: str,
        variants: List[Variant],
        bam: AlignmentFile,
        min_base_quality: int,
        min_mapping_quality: int,
        include_ambiguous_bases: bool = False) -> Dict[Tuple, 'CoverageMetrics']:
    """
    Opens ONE pileup iterator over the entire chromosome region covered by variants.
    Metrics are computed IMMEDIATELY for each pileup column while it is still valid —
    avoids segfaults from storing PileupColumn objects after the iterator advances.

    Returns {(pos, REF, ALT[0]): CoverageMetrics}.
    """
    if not variants:
        return {}

    # index variants by 1-based position for O(1) lookup during streaming
    variants_by_pos: Dict[int, List[Variant]] = defaultdict(list)
    for v in variants:
        variants_by_pos[v.POS].append(v)

    start = variants[0].POS - 1   # 0-based inclusive
    end = variants[-1].POS        # exclusive end for pysam

    results: Dict[Tuple, CoverageMetrics] = {}

    for pileup_col in get_region_pileup(
            chrom=chrom, start=start, end=end,
            bam=bam,
            min_base_quality=min_base_quality,
            min_mapping_quality=min_mapping_quality,
    ):
        ref_pos = pileup_col.reference_pos + 1  # convert to 1-based
        if ref_pos not in variants_by_pos:
            continue

        # compute metrics NOW while pileup_col is still valid in C memory
        for variant in variants_by_pos[ref_pos]:
            metrics = _get_metrics_from_column(variant, pileup_col, include_ambiguous_bases)
            if metrics is not None:
                results[(ref_pos, variant.REF, variant.ALT[0])] = metrics

    return results


@dataclass
class CoverageMetrics:
    # number supporting reads of each base, including the reference
    ac: dict
    # total depth of coverage
    dp: int
    # median base call quality of each base, including the reference
    bqs: dict = None
    # median mapping quality of each alternate base, including the reference
    mqs: dict = None
    # median position within the read of each alternate base, including the reference
    positions: dict = None
    # base call quality distribution of each base, including the reference
    all_bqs: dict = None
    # mapping quality distribution of each base, including the reference
    all_mqs: dict = None
    # position within the read distribution of each base, including the reference
    all_positions: dict = None


EMPTY_METRICS = CoverageMetrics(
    ac=Counter(), dp=0, bqs=Counter(), mqs=Counter(), positions=Counter(),
    all_bqs={}, all_mqs={}, all_positions={}
)


def _get_metrics_from_column(variant: Variant, pileup_col,
                              include_ambiguous_bases=False) -> 'CoverageMetrics':
    """Dispatch to the right metrics function based on variant type."""
    if is_snp(variant):
        return _get_snv_metrics_from_column(pileup_col, include_ambiguous_bases)
    elif is_insertion(variant):
        return _get_insertion_metrics_from_column(variant, pileup_col)
    elif is_deletion(variant):
        return _get_deletion_metrics_from_column(variant, pileup_col)
    return None


def _get_snv_metrics_from_column(pileup_col, include_ambiguous_bases=False) -> CoverageMetrics:
    bases = []
    qualities = []
    mapping_qualities = []
    query_positions = []

    for read in pileup_col.pileups:
        if read.is_refskip:
            continue
        if read.is_del:
            bases.append("")
            qualities.append(0)
            mapping_qualities.append(read.alignment.mapping_quality)
            query_positions.append(read.query_position_or_next)
        else:
            base = read.alignment.query_sequence[read.query_position].upper()
            bases.append(base)
            qualities.append(read.alignment.query_qualities[read.query_position])
            mapping_qualities.append(read.alignment.mapping_quality)
            query_positions.append(read.query_position)

    all_bqs = aggregate_list_per_base(bases, qualities)
    all_mqs = aggregate_list_per_base(bases, mapping_qualities)
    all_positions = aggregate_list_per_base(bases, query_positions)

    bqs = Counter({b: np.median(l) for b, l in all_bqs.items()})
    mqs = Counter({b: np.median(l) for b, l in all_mqs.items()})
    positions = Counter({b: np.median(l) for b, l in all_positions.items()})

    ac = Counter(b for b in bases if b != "")

    if include_ambiguous_bases:
        dp = len(bases)
    else:
        dp = sum(1 for b in bases if b == "" or b not in AMBIGUOUS_BASES)

    return CoverageMetrics(
        ac=ac, dp=dp, bqs=bqs, mqs=mqs, positions=positions,
        all_bqs=all_bqs, all_mqs=all_mqs, all_positions=all_positions
    )


def _get_insertion_metrics_from_column(variant: Variant, pileup_col) -> CoverageMetrics:
    ac = {alt.upper(): 0 for alt in variant.ALT}
    mq = {alt.upper(): [] for alt in variant.ALT}
    mq[variant.REF] = []
    pos = {alt.upper(): [] for alt in variant.ALT}
    pos[variant.REF] = []
    variant_position = variant.POS
    insertion_length = len(variant.ALT[0]) - len(variant.REF)
    insertion = variant.ALT[0][1:]
    alt_upper = variant.ALT[0].upper()

    pileup_reads = pileup_col.pileups
    dp = len(pileup_reads)

    for pileup_read in pileup_reads:
        if pileup_read.indel > 0:
            index = pileup_read.alignment.reference_start
            relative_position = 0
            for cigar_type, cigar_length in pileup_read.alignment.cigartuples:
                if cigar_type in [0, 2, 3, 7, 8]:
                    index += cigar_length
                    if index > variant_position:
                        break
                if cigar_type in [0, 1, 4, 7, 8]:
                    relative_position += cigar_length
                if cigar_type == 1:
                    insertion_in_query = pileup_read.alignment.query[
                                         relative_position:relative_position + insertion_length]
                    if index == variant_position and cigar_length == insertion_length \
                            and insertion == insertion_in_query:
                        ac[alt_upper] += 1
                        mq[alt_upper].append(pileup_read.alignment.mapping_quality)
                        pos[alt_upper].append(pileup_read.query_position_or_next)
        elif pileup_read.indel == 0:
            mq[variant.REF].append(pileup_read.alignment.mapping_quality)
            pos[variant.REF].append(pileup_read.query_position_or_next)

    return CoverageMetrics(
        ac=Counter(ac), dp=dp,
        mqs=Counter({k: np.median(l) for k, l in mq.items()}),
        positions=Counter({k: np.median(l) for k, l in pos.items()}),
        bqs=Counter(),
        all_mqs={k: l for k, l in mq.items()},
        all_positions={k: l for k, l in pos.items()},
        all_bqs=Counter()
    )


def _get_deletion_metrics_from_column(variant: Variant, pileup_col) -> CoverageMetrics:
    ac = {alt.upper(): 0 for alt in variant.ALT}
    mq = {alt.upper(): [] for alt in variant.ALT}
    mq[variant.REF] = []
    pos = {alt.upper(): [] for alt in variant.ALT}
    pos[variant.REF] = []
    variant_position = variant.POS
    deletion_length = len(variant.REF) - len(variant.ALT[0])
    alt_upper = variant.ALT[0].upper()

    pileup_reads = pileup_col.pileups
    dp = len(pileup_reads)

    for pileup_read in pileup_reads:
        if pileup_read.indel < 0:
            start = pileup_read.alignment.reference_start
            for cigar_type, cigar_length in pileup_read.alignment.cigartuples:
                if cigar_type in [0, 3, 7, 8]:
                    start += cigar_length
                elif cigar_type == 2:
                    if start == variant_position and cigar_length == deletion_length:
                        ac[alt_upper] += 1
                        mq[alt_upper].append(pileup_read.alignment.mapping_quality)
                        pos[alt_upper].append(pileup_read.query_position_or_next)
                        break
                    else:
                        start += cigar_length
                if start > variant_position:
                    break
        elif pileup_read.indel == 0:
            mq[variant.REF].append(pileup_read.alignment.mapping_quality)
            pos[variant.REF].append(pileup_read.query_position_or_next)

    return CoverageMetrics(
        ac=Counter(ac), dp=dp,
        mqs=Counter({k: np.median(l) for k, l in mq.items()}),
        positions=Counter({k: np.median(l) for k, l in pos.items()}),
        bqs=Counter(),
        all_mqs={k: l for k, l in mq.items()},
        all_positions={k: l for k, l in pos.items()},
        all_bqs=Counter()
    )

def aggregate_list_per_base(bases, values) -> dict:
    aggregated_values = {}
    for b, v in zip(bases, values):
        if b not in aggregated_values:
            aggregated_values[b] = []
        aggregated_values[b].append(v)
    return aggregated_values