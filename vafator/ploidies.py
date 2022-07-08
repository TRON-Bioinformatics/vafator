import os
from typing import Union

from cyvcf2 import Variant
from pybedtools import BedTool, Interval

from vafator.tests.utils import VafatorVariant

DEFAULT_PLOIDY = 2.0


class PloidyManager:

    def __init__(self, local_copy_numbers: str = None, genome_wide_ploidy: float = DEFAULT_PLOIDY):

        if local_copy_numbers is not None and not os.path.exists(local_copy_numbers):
            raise ValueError('The provided tumor ploidy is neither a copy number value or a BED file with copy '
                             'numbers')
        self.bed = BedTool(local_copy_numbers) if local_copy_numbers is not None else None
        self.ploidy = genome_wide_ploidy

    def get_ploidy(self, variant: Union[Variant, VafatorVariant]) -> float:

        result = self.ploidy
        if self.bed is not None:
            # read from the BED file
            # NOTE: converts from 1-based into 0-based
            hits = self.bed.all_hits(Interval(chrom=variant.CHROM, start=variant.POS - 1, end=variant.POS - 1))
            if len(hits) > 0:
                result = float(hits[0].score)

        return result


default_ploidy_manager = PloidyManager()
