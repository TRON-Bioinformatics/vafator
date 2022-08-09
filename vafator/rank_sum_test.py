from typing import List, Tuple
import scipy.stats
import numpy as np


def calculate_rank_sum_test(alternate_dist: List[int], reference_dist: List[int]) -> Tuple[float, float]:
    stat, pvalue = scipy.stats.ranksums(x=alternate_dist, y=reference_dist)
    return round(stat, 3), round(pvalue, 5)


def get_rank_sum_tests(distributions: dict, variant):
    stats = []
    pvalues = []
    for alt in variant.ALT:
        stat, pvalue = calculate_rank_sum_test(
            alternate_dist=distributions.get(alt, []),
            reference_dist=distributions.get(variant.REF, []))
        if not np.isnan(stat) and not np.isnan(pvalue):
            stats.append(str(stat))
            pvalues.append(str(pvalue))
    return pvalues, stats


