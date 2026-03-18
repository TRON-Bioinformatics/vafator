from unittest import TestCase
import pkg_resources

from vafator.tests.utils import VafatorVariant
from vafator.ploidies import PloidyManager, default_ploidy_manager


class PloidyManagerTest(TestCase):

    def test_default_ploidy_manager(self):
        self.assertEqual(PloidyManager().get_ploidy(
            variant=VafatorVariant(chromosome="chr1", position=12345, reference="A", alternative="C")), 2.0)
        self.assertEqual(default_ploidy_manager.get_ploidy(
            variant=VafatorVariant(chromosome="chr1", position=12345, reference="A", alternative="C")), 2.0)
        self.assertEqual(default_ploidy_manager.get_ploidy(variant=None), 2.0)

    def test_genome_wide_ploidy_manager(self):
        self.assertEqual(PloidyManager(genome_wide_ploidy=3.2).get_ploidy(
            variant=VafatorVariant(chromosome="chr1", position=12345, reference="A", alternative="C")), 3.2)
        self.assertEqual(PloidyManager(genome_wide_ploidy=3.2).get_ploidy(variant=None), 3.2)

    def test_local_copy_numbers_ploidy_manager(self):
        input_bed = pkg_resources.resource_filename(__name__, "resources/test_copy_numbers.bed")
        self.assertEqual(PloidyManager(local_copy_numbers=input_bed).get_ploidy(
            variant=VafatorVariant(chromosome="chr1", position=12345, reference="A", alternative="C")), 1.2)
        self.assertEqual(PloidyManager(local_copy_numbers=input_bed).get_ploidy(
            variant=VafatorVariant(chromosome="chr2", position=12345, reference="A", alternative="C")), 3.2)
        # test non existing interval
        self.assertEqual(PloidyManager(local_copy_numbers=input_bed).get_ploidy(
            variant=VafatorVariant(chromosome="chr3", position=12345, reference="A", alternative="C")), 2.0)

    def test_interval_boundaries(self):
        input_bed = pkg_resources.resource_filename(__name__, "resources/test_copy_numbers.bed")
        # lower boundary — POS 10000 is 0-based 9999, outside interval start 10000
        self.assertEqual(PloidyManager(local_copy_numbers=input_bed).get_ploidy(
            variant=VafatorVariant(chromosome="chr1", position=10000, reference="A", alternative="C")), 2.0)
        # just inside lower boundary
        self.assertEqual(PloidyManager(local_copy_numbers=input_bed).get_ploidy(
            variant=VafatorVariant(chromosome="chr1", position=10001, reference="A", alternative="C")), 1.2)
        # upper boundary
        self.assertEqual(PloidyManager(local_copy_numbers=input_bed).get_ploidy(
            variant=VafatorVariant(chromosome="chr1", position=20000, reference="A", alternative="C")), 1.2)
        self.assertEqual(PloidyManager(local_copy_numbers=input_bed).get_ploidy(
            variant=VafatorVariant(chromosome="chr1", position=20001, reference="A", alternative="C")), 2.1)

    def test_invalid_bed_raises(self):
        with self.assertRaises(ValueError):
            PloidyManager(local_copy_numbers="/nonexistent/path.bed")
