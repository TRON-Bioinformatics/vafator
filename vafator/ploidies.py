import os
from typing import Union
import pandas as pd
from cyvcf2 import Variant

from vafator.tests.utils import VafatorVariant

DEFAULT_PLOIDY = 2.0


class PloidyManager:

    def __init__(self, local_copy_numbers: str = None, genome_wide_ploidy: float = DEFAULT_PLOIDY):

        if local_copy_numbers is not None and not os.path.exists(local_copy_numbers):
            raise ValueError('The provided tumor ploidy is neither a copy number value or a BED file with copy '
                             'numbers')
        self.report_value = local_copy_numbers if local_copy_numbers else genome_wide_ploidy
        self.bed = pd.read_csv(local_copy_numbers, sep='\t', names=['chromosome', 'start', 'end', 'copy_number']) \
            if local_copy_numbers is not None else None
        self.ploidy = genome_wide_ploidy

    def get_ploidy(self, variant: Union[Variant, VafatorVariant]) -> float:

        result = self.ploidy
        if self.bed is not None:
            # read from the BED file
            # NOTE: converts variant position from 1-based into 0-based and considers intervals as half-closed
            hits = self.bed[(self.bed.chromosome == variant.CHROM) &
                            (self.bed.start <= variant.POS - 1) &
                            (self.bed.end > variant.POS - 1)]
            if hits.shape[0] > 0:
                result = float(hits.copy_number.iloc[0])

        return result


default_ploidy_manager = PloidyManager()
