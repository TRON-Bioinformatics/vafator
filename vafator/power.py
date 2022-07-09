from scipy.stats import binom
from vafator.ploidies import default_ploidy_manager


DEFAULT_PURITY = 1.0
DEFAULT_NORMAL_PLOIDY = 2


class PowerCalculator:

    def __init__(self, tumor_ploidies, purities, normal_ploidy=DEFAULT_NORMAL_PLOIDY):
        self.normal_ploidy = normal_ploidy
        self.purities = purities
        self.tumor_ploidies = tumor_ploidies

    def calculate_power(self, dp, ac, sample, variant):
        """
        Return the binomial probability of observing ac supporting reads, given a total coverage dp and a
        expected VAF tumor purity / 2.
        """
        # NOTE: assumes normal ploidy of 2, this will not hold in sexual chromosomes except PARs or other no diploid
        # organisms
        expected_vaf = self.calculate_expected_vaf(sample, variant)
        pvalue = binom.cdf(k=ac, n=dp, p=expected_vaf)
        return pvalue

    def calculate_expected_vaf(self, sample, variant):
        purity = self.purities.get(sample, DEFAULT_PURITY)
        tumor_ploidy = self.tumor_ploidies.get(sample, default_ploidy_manager).get_ploidy(variant=variant)
        corrected_tumor_ploidy = purity * tumor_ploidy + ((1 - purity) * self.normal_ploidy)
        expected_vaf = purity / corrected_tumor_ploidy
        return expected_vaf
