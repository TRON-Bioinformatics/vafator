from typing import List, Tuple
import scipy.stats


def calculate_rank_sum_test(alternate_dist: List[int], reference_dist: List[int]) -> Tuple[float, float]:
    stat, pvalue = scipy.stats.ranksums(x=alternate_dist, y=reference_dist)
    return round(stat, 3), round(pvalue, 5)


