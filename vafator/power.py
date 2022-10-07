from typing import Optional

from cyvcf2 import Variant
from scipy.stats import binom
from vafator.ploidies import default_ploidy_manager


DEFAULT_PURITY = 1.0
DEFAULT_NORMAL_PLOIDY = 2
DEFAULT_FPR = 5*(10**-7)
DEFAULT_ERROR_RATE = 10**-3


class PowerCalculator:

    def __init__(
            self,
            tumor_ploidies: dict,
            purities: dict,
            normal_ploidy: int = DEFAULT_NORMAL_PLOIDY,
            fpr: float = DEFAULT_FPR,
            error_rate: float = DEFAULT_ERROR_RATE):

        self.normal_ploidy = normal_ploidy
        self.purities = purities
        self.tumor_ploidies = tumor_ploidies
        self.fpr = fpr
        self.error_rate = error_rate

    def calculate_power(self, dp: int, ac: int, sample: str, variant: Optional[Variant]) -> float:
        """
        Return the binomial probability of observing ac or less supporting reads, given a total coverage dp and a
        expected VAF tumor purity / 2.
        """
        # NOTE: assumes normal ploidy of 2, this will not hold in sexual chromosomes except PARs or other no diploid
        # organisms
        expected_vaf = self.calculate_expected_vaf(sample, variant)
        pvalue = binom.cdf(k=ac, n=dp, p=expected_vaf)
        return round(pvalue, 5)

    def calculate_expected_vaf(self, sample: str, variant: Optional[Variant]) -> float:
        """
        Calculates expected VAF based on purity and local copy number.
        expected VAF = purity / adjusted tumor ploidy

        In a scenario with purity = 1, tumor CN = 2 and normal CN = 2 => expected VAF = 0.5
        """
        purity = self.purities.get(sample, DEFAULT_PURITY)
        tumor_ploidy = max(1, self.tumor_ploidies.get(sample, default_ploidy_manager).get_ploidy(variant=variant))
        corrected_tumor_ploidy = purity * tumor_ploidy + ((1 - purity) * self.normal_ploidy)
        expected_vaf = purity / corrected_tumor_ploidy
        return expected_vaf

    def _calculate_p(self, m: int, n: int) -> float:
        """
        This is the P defined in Carter, 2012, where m is the number of observed reads supporting the mutation and
        n is the total coverage (dp).
        This calculates the probability of observing k or more identical non-reference reads due to sequencing error.
        """
        result = 1
        if m >= 1:
            result = 1 - binom.cdf(k=m - 1, n=n, p=self.error_rate / 3)
        return result

    def _calculate_d(self, k: int, n: int) -> float:
        """
        This is the d defined in Carter, 2012, where k is the number of observed reads supporting the mutation and
        n is the total coverage (dp)
        Ratio between difference between FPR and prob. of observing >= k reads and
        difference between prob. or observing k - 1 reads and observing k reads
        """
        p = self._calculate_p(m=k, n=n)
        p_1 = self._calculate_p(m=k - 1, n=n)
        return (self.fpr - p) / (p_1 - p)

    def _calculate_k(self, dp: int) -> int:
        k = 1
        while self._calculate_p(m=k, n=dp) > self.fpr:
            k += 1

        return k

    def calculate_absolute_power(self, sample, variant, dp: int) -> float:
        """
        This is the power as defined in Carter, 2012, where n is the total coverage (dp),
        f is the expected variant allele frequency and k is minimum number of reads supporting the mutation where
        P <= FPR
        """
        k = self._calculate_k(dp=dp)
        n = dp
        f = self.calculate_expected_vaf(sample, variant)
        power = 1 - binom.cdf(k=k - 1, n=n, p=f) + self._calculate_d(k=k, n=n) * binom(n, f).pmf(k)
        return round(power, 5), k


