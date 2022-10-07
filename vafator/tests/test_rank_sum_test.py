from unittest import TestCase
from vafator.rank_sum_test import calculate_rank_sum_test


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
