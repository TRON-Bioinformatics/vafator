from math import isnan
from unittest import TestCase
from vafator.rank_sum_test import calculate_rank_sum_test, get_rank_sum_tests
from vafator.pileups import VariantRecord


class TestRankSumTest(TestCase):

    def test_equal_distribution(self):
        distribution = [1, 2, 3, 4, 5, 6, 7, 8, 9]
        stat, pvalue = calculate_rank_sum_test(distribution, distribution)
        self.assertEqual(stat, 0.0)
        self.assertEqual(pvalue, 1.0)

    def test_direction(self):
        distribution1 = [1, 2, 3, 4, 5, 6]
        distribution2 = [4, 5, 6, 7, 8, 9]
        stat1, pvalue1 = calculate_rank_sum_test(distribution1, distribution2)
        self.assertLess(stat1, 0.0)
        self.assertLess(pvalue1, 1.0)
        stat2, pvalue2 = calculate_rank_sum_test(distribution2, distribution1)
        self.assertGreater(stat2, 0.0)
        self.assertLess(pvalue2, 1.0)

        self.assertEqual(stat1, -stat2)

    def test_gatk_example(self):
        """
        This test data comes from https://gatk.broadinstitute.org/hc/en-us/articles/360035531952-Rank-Sum-Test
        """
        ref_distribution = [20, 25, 26, 30, 32, 40, 47, 50, 53, 60]
        alt_distribution = [0, 7, 10, 17, 20, 21, 30, 34, 40, 45]
        stat, pvalue = calculate_rank_sum_test(alt_distribution, ref_distribution)
        self.assertEqual(stat, -2.154)
        self.assertEqual(pvalue, 0.03121)

    def test_empty_alternate_returns_nan(self):
        stat, pvalue = calculate_rank_sum_test([], [1, 2, 3])
        self.assertTrue(isnan(stat))
        self.assertTrue(isnan(pvalue))

    def test_empty_reference_returns_nan(self):
        stat, pvalue = calculate_rank_sum_test([1, 2, 3], [])
        self.assertTrue(isnan(stat))
        self.assertTrue(isnan(pvalue))

    def test_both_empty_returns_nan(self):
        stat, pvalue = calculate_rank_sum_test([], [])
        self.assertTrue(isnan(stat))
        self.assertTrue(isnan(pvalue))

    def test_get_rank_sum_tests_snv(self):
        variant = VariantRecord(CHROM='chr1', POS=100, REF='A', ALT=['T'])
        distributions = {
            'A': [20, 25, 30, 35, 40],  # ref is higher
            'T': [1, 5, 10, 15, 20],    # alt is lower
        }
        pvalues, stats = get_rank_sum_tests(distributions, variant)
        self.assertEqual(len(stats), 1)
        self.assertEqual(len(pvalues), 1)
        # alt < ref so ranksums(alt, ref) returns negative stat
        self.assertLess(float(stats[0]), 0.0)

    def test_get_rank_sum_tests_no_alt_reads(self):
        variant = VariantRecord(CHROM='chr1', POS=100, REF='A', ALT=['T'])
        distributions = {'A': [20, 25, 30]}  # no T reads
        pvalues, stats = get_rank_sum_tests(distributions, variant)
        self.assertEqual(stats, [])
        self.assertEqual(pvalues, [])