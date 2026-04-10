from collections import Counter
from dataclasses import dataclass
from math import nan
from typing import List, Dict

import numpy as np


@dataclass
class VariantRecord:
    """Lightweight, picklable variant representation used by pileup workers.
    Mirrors the cyvcf2.Variant fields accessed by pileup and metrics functions.

    Attributes:
        CHROM: chromosome name
        POS: 1-based variant position
        REF: reference allele
        ALT: list of alternate alleles
    """

    CHROM: str
    POS: int
    REF: str
    ALT: List[str]


@dataclass
class CoverageMetrics:
    """Pileup metrics computed for a single variant in a single BAM.

    Attributes:
        ac: allele counts per base, including the reference
        dp: total depth of coverage
        bqs: median base call quality per allele, including the reference
        mqs: median mapping quality per allele, including the reference
        positions: median read position per allele, including the reference
        all_bqs: full base call quality distribution per allele
        all_mqs: full mapping quality distribution per allele
        all_positions: full read position distribution per allele
    """

    ac: dict
    dp: int
    bqs: dict = None
    mqs: dict = None
    positions: dict = None
    all_bqs: dict = None
    all_mqs: dict = None
    all_positions: dict = None


EMPTY_METRICS = CoverageMetrics(
    ac=Counter(),
    dp=0,
    bqs=Counter(),
    mqs=Counter(),
    positions=Counter(),
    all_bqs={},
    all_mqs={},
    all_positions={},
)


def aggregate_list_per_base(bases: List[str], values: list) -> Dict[str, list]:
    """Group a list of values by their corresponding base.

    Args:
        bases: list of base characters (e.g. ['A', 'T', 'A', ''])
        values: list of numeric values parallel to bases

    Returns:
        dict mapping each base to a list of its associated values
    """
    aggregated_values = {}
    for b, v in zip(bases, values):
        if b not in aggregated_values:
            aggregated_values[b] = []
        aggregated_values[b].append(v)
    return aggregated_values


def safe_median(values: list) -> float:
    """Return the median of a list, or nan if the list is empty.
    Avoids numpy RuntimeWarning raised by np.median on empty arrays.

    Args:
        values: list of numeric values

    Returns:
        median as float, or nan if values is empty
    """
    return float(np.median(values)) if values else nan


def is_snp(variant) -> bool:
    """Return True if the variant is a single nucleotide polymorphism.

    Args:
        variant: any object with REF and ALT attributes (Variant or VariantRecord)
    """
    return len(variant.REF) == 1 and len(variant.ALT[0]) == 1


def is_insertion(variant) -> bool:
    """Return True if the variant is an insertion.

    Args:
        variant: any object with REF and ALT attributes (Variant or VariantRecord)
    """
    return len(variant.REF) == 1 and len(variant.ALT[0]) > 1


def is_deletion(variant) -> bool:
    """Return True if the variant is a deletion.

    Args:
        variant: any object with REF and ALT attributes (Variant or VariantRecord)
    """
    return len(variant.ALT[0]) == 1 and len(variant.REF) > 1
