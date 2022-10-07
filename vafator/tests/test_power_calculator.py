from unittest import TestCase
import pkg_resources

from vafator.ploidies import PloidyManager
from vafator.power import PowerCalculator
from vafator.tests.utils import VafatorVariant


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

    def test_eaf_copy_number_below_one(self):
        power = PowerCalculator(
            tumor_ploidies={'tumor': PloidyManager(genome_wide_ploidy=0.5)}, purities={'tumor': 0.9})
        self.assertLessEqual(power.calculate_expected_vaf(sample='tumor', variant=None), 1.0)

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
        ploidy_manager = {'tumor': PloidyManager(genome_wide_ploidy=2)}

        calculator = PowerCalculator(tumor_ploidies=ploidy_manager, purities={'tumor': 0.8})
        p, k = calculator.calculate_absolute_power(dp=100, sample='tumor', variant=None)
        self.assertEqual(p, 1.0)
        self.assertEqual(k, 4)
        p, k = calculator.calculate_absolute_power(dp=50, sample='tumor', variant=None)
        self.assertEqual(p, 1.0)
        self.assertEqual(k, 4)
        p, k = calculator.calculate_absolute_power(dp=10, sample='tumor', variant=None)
        self.assertEqual(p, 0.85408)
        self.assertEqual(k, 3)
        p, k = calculator.calculate_absolute_power(dp=2, sample='tumor', variant=None)
        self.assertEqual(p, 0.16009)
        self.assertEqual(k, 2)

        calculator = PowerCalculator(tumor_ploidies=ploidy_manager, purities={'tumor': 0.6})
        p, k = calculator.calculate_absolute_power(
            dp=100, sample='tumor', variant=None)
        self.assertEqual(p, 1.0)
        self.assertEqual(k, 4)
        p, k = calculator.calculate_absolute_power(dp=50, sample='tumor', variant=None)
        self.assertEqual(p, 1.00007)
        self.assertEqual(k, 4)
        p, k = calculator.calculate_absolute_power(dp=10, sample='tumor', variant=None)
        self.assertEqual(p, 0.64373)
        self.assertEqual(k, 3)
        p, k = calculator.calculate_absolute_power(dp=2, sample='tumor', variant=None)
        self.assertEqual(p, 0.09005)
        self.assertEqual(k, 2)

        calculator = PowerCalculator(tumor_ploidies=ploidy_manager, purities={'tumor': 0.1})
        p, k = calculator.calculate_absolute_power(
            dp=100, sample='tumor', variant=None)
        self.assertEqual(p, 0.75607)
        self.assertEqual(k, 4)
        p, k = calculator.calculate_absolute_power(dp=50, sample='tumor', variant=None)
        self.assertEqual(p, 0.33419)
        self.assertEqual(k, 4)
        p, k = calculator.calculate_absolute_power(dp=10, sample='tumor', variant=None)
        self.assertEqual(p, 0.01254)
        self.assertEqual(k, 3)
        p, k = calculator.calculate_absolute_power(dp=2, sample='tumor', variant=None)
        self.assertEqual(p, 0.0025)
        self.assertEqual(k, 2)
