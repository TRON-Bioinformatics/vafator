from unittest import TestCase

import numpy as np
import pkg_resources

from ploidies import PloidyManager
from power import PowerCalculator
from tests.utils import VafatorVariant


class PowerCalculatorTest(TestCase):

    def test_power_calculator(self):
        power = PowerCalculator(
            tumor_ploidies={'tumor': PloidyManager(genome_wide_ploidy=2.5)}, purities={'tumor': 0.8})
        self.assertAlmostEqual(power.calculate_power(dp=10, ac=0, sample='tumor', variant=None), 0.01734)
        self.assertAlmostEqual(power.calculate_power(dp=10, ac=1, sample='tumor', variant=None), 0.10404917949499565, 5)
        self.assertAlmostEqual(power.calculate_power(dp=10, ac=2, sample='tumor', variant=None), 0.29914139104811255, 5)
        self.assertAlmostEqual(power.calculate_power(dp=10, ac=3, sample='tumor', variant=None), 0.5592643397856016, 5)
        self.assertAlmostEqual(power.calculate_power(dp=10, ac=4, sample='tumor', variant=None), 0.7868719199309048, 5)
        self.assertAlmostEqual(power.calculate_power(dp=10, ac=5, sample='tumor', variant=None), 0.9234364680180867, 5)
        self.assertAlmostEqual(power.calculate_power(dp=10, ac=6, sample='tumor', variant=None), 0.9803383630544125, 5)
        self.assertAlmostEqual(power.calculate_power(dp=10, ac=7, sample='tumor', variant=None), 0.9965960473505056, 5)
        self.assertAlmostEqual(power.calculate_power(dp=10, ac=8, sample='tumor', variant=None), 0.999644363156023, 5)
        self.assertAlmostEqual(power.calculate_power(dp=10, ac=9, sample='tumor', variant=None), 0.9999830649121916, 5)
        self.assertAlmostEqual(power.calculate_power(dp=10, ac=10, sample='tumor', variant=None), 1.0)
        self.assertAlmostEqual(power.calculate_power(dp=10, ac=11, sample='tumor', variant=None), 1.0)

    def test_varying_purity(self):
        power1 = PowerCalculator(
            tumor_ploidies={'tumor': PloidyManager(genome_wide_ploidy=2.5)}, purities={'tumor': 0.8})
        power2 = PowerCalculator(
            tumor_ploidies={'tumor': PloidyManager(genome_wide_ploidy=2.5)}, purities={'tumor': 0.6})
        self.assertLess(
            power1.calculate_power(dp=10, ac=2, sample='tumor', variant=None),
            power2.calculate_power(dp=10, ac=2, sample='tumor', variant=None))

        power3 = PowerCalculator(
            tumor_ploidies={'tumor': PloidyManager(genome_wide_ploidy=2.5)}, purities={'tumor': 0.4})
        self.assertLess(
            power2.calculate_power(dp=10, ac=2, sample='tumor', variant=None),
            power3.calculate_power(dp=10, ac=2, sample='tumor', variant=None))

    def test_varying_ploidy(self):
        power1 = PowerCalculator(
            tumor_ploidies={'tumor': PloidyManager(genome_wide_ploidy=2.0)}, purities={'tumor': 0.8})
        power2 = PowerCalculator(
            tumor_ploidies={'tumor': PloidyManager(genome_wide_ploidy=4.0)}, purities={'tumor': 0.8})
        self.assertLess(
            power1.calculate_power(dp=10, ac=2, sample='tumor', variant=None),
            power2.calculate_power(dp=10, ac=2, sample='tumor', variant=None))

        power3 = PowerCalculator(
            tumor_ploidies={'tumor': PloidyManager(genome_wide_ploidy=6.0)}, purities={'tumor': 0.8})
        self.assertLess(
            power2.calculate_power(dp=10, ac=2, sample='tumor', variant=None),
            power3.calculate_power(dp=10, ac=2, sample='tumor', variant=None))

    def test_local_copy_numbers(self):
        input_bed = pkg_resources.resource_filename(__name__, "resources/test_copy_numbers.bed")
        power = PowerCalculator(
            tumor_ploidies={'tumor': PloidyManager(local_copy_numbers=input_bed)}, purities={'tumor': 0.8})

        self.assertNotEqual(
            power.calculate_power(dp=10, ac=2, sample='tumor',
                                  variant=VafatorVariant(
                                      chromosome='chr1', position=10001, reference='A', alternative='G')),
            power.calculate_power(dp=10, ac=2, sample='tumor',
                                  variant=VafatorVariant(
                                      chromosome='chr1', position=20001, reference='A', alternative='G')))
        self.assertEqual(
            power.calculate_power(dp=10, ac=2, sample='tumor',
                                  variant=VafatorVariant(
                                      chromosome='chr1', position=10001, reference='A', alternative='G')),
            power.calculate_power(dp=10, ac=2, sample='tumor',
                                  variant=VafatorVariant(
                                      chromosome='chr1', position=10002, reference='A', alternative='G')))

    def test_absolute_power_calculator(self):
        power = PowerCalculator(
            tumor_ploidies={'tumor': PloidyManager(genome_wide_ploidy=2.5)}, purities={'tumor': 0.8})
        self.assertAlmostEqual(power.calculate_absolute_power(dp=10, ac=0, sample='tumor', variant=None), 0.0, 5)
        self.assertAlmostEqual(power.calculate_absolute_power(dp=10, ac=1, sample='tumor', variant=None), 0.9823689574968194, 5)
        self.assertAlmostEqual(power.calculate_absolute_power(dp=10, ac=2, sample='tumor', variant=None), 0.8956871759303671, 5)
        self.assertAlmostEqual(power.calculate_absolute_power(dp=10, ac=3, sample='tumor', variant=None), 0.7267089415880859, 5)
        self.assertAlmostEqual(power.calculate_absolute_power(dp=10, ac=4, sample='tumor', variant=None), 26.106282321427628, 5)
        self.assertAlmostEqual(power.calculate_absolute_power(dp=10, ac=5, sample='tumor', variant=None), 26389.64274484705, 5)
        self.assertAlmostEqual(power.calculate_absolute_power(dp=10, ac=6, sample='tumor', variant=None), 28473705.99681597, 5)
        self.assertAlmostEqual(power.calculate_absolute_power(dp=10, ac=7, sample='tumor', variant=None), np.inf)
        self.assertAlmostEqual(power.calculate_absolute_power(dp=10, ac=8, sample='tumor', variant=None), np.inf)
        self.assertAlmostEqual(power.calculate_absolute_power(dp=10, ac=9, sample='tumor', variant=None), np.inf)
        self.assertAlmostEqual(power.calculate_absolute_power(dp=10, ac=10, sample='tumor', variant=None), np.inf)
        self.assertTrue(np.isnan(power.calculate_absolute_power(dp=10, ac=11, sample='tumor', variant=None)))