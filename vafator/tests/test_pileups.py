from collections import Counter
from unittest import TestCase
import pkg_resources
import pysam

from vafator.tests.utils import VafatorVariant
from vafator.pileups import get_variant_pileup, get_metrics


class TestPileups(TestCase):

    min_base_quality = 0
    min_mapping_quality = 0
    bam_file = pkg_resources.resource_filename(
        __name__,
        "resources/project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.chr1_1000000_2000000.bam")
    bam_reader = pysam.AlignmentFile(bam_file)

    def test_snv_metrics(self):
        variant = VafatorVariant(chromosome="chr1", position=1017341, reference="G", alternative=["T"])
        self._assert_metrics(variant=variant, expected_ac={'G': 5, 'T': 6}, expected_dp=11)
        self._assert_metrics(variant=variant, expected_ac={'G': 1, 'T': 2}, expected_dp=3, min_base_quality=40)

    def test_snv_metrics_2(self):
        variant = VafatorVariant(chromosome="chr1", position=1018144, reference="T", alternative=["C"])
        self._assert_metrics(variant=variant, expected_ac={'C': 9, 'T': 11}, expected_dp=20)
        self._assert_metrics(variant=variant, expected_ac={'C': 3, 'T': 4}, expected_dp=7, min_base_quality=40)
        self._assert_metrics(variant=variant, expected_ac=Counter(), expected_dp=0, min_mapping_quality=65)

    def test_insertion_metrics(self):
        # variant called in the VCF shows no read support (!?)
        variant = VafatorVariant(chromosome="chr1", position=1247578, reference="T", alternative=["TGG"])
        self._assert_metrics(variant=variant, expected_ac={'TGG': 0}, expected_dp=3)
        # there is one read supporting this insertion of 3 Gs
        variant = VafatorVariant(chromosome="chr1", position=1247578, reference="T", alternative=["TGGG"])
        self._assert_metrics(variant=variant, expected_ac={'TGGG': 1}, expected_dp=3)
        # this ensures that the insertion sequence is checked not only the insertion length!
        variant = VafatorVariant(chromosome="chr1", position=1247578, reference="T", alternative=["TGGA"])
        self._assert_metrics(variant=variant, expected_ac={'TGGA': 0}, expected_dp=3)
        # there is one read supporting this insertion of 4 Gs
        variant = VafatorVariant(chromosome="chr1", position=1247578, reference="T", alternative=["TGGGG"])
        self._assert_metrics(variant=variant, expected_ac={'TGGGG': 1}, expected_dp=3)
        # there is no read supporting an insertion of 5 Gs
        variant = VafatorVariant(chromosome="chr1", position=1247578, reference="T", alternative=["TGGGGG"])
        self._assert_metrics(variant=variant, expected_ac={'TGGGGG': 0}, expected_dp=3)

    def test_insertion_metrics_2(self):
        variant = VafatorVariant(chromosome="chr1", position=1594199, reference="C", alternative=["CT"])
        self._assert_metrics(variant=variant, expected_ac={'CT': 9}, expected_dp=11)
        self._assert_metrics(variant=variant, expected_ac={'CT': 5}, expected_dp=7, min_mapping_quality=40)
        self._assert_metrics(variant=variant, expected_ac={'CT': 4}, expected_dp=4, min_base_quality=40)

    def test_deletion_metrics(self):
        variant = VafatorVariant(chromosome="chr1", position=1510035, reference="GGC", alternative=["G"])
        self._assert_metrics(variant=variant, expected_ac={'G': 12}, expected_dp=13)
        self._assert_metrics(variant=variant, expected_ac={'G': 0}, expected_dp=0, min_mapping_quality=61)
        self._assert_metrics(variant=variant, expected_ac={'G': 3}, expected_dp=3, min_base_quality=40)
        # deletions with a reference sequence not matching the reference would be matched
        # vafator expects correct indel calls
        variant = VafatorVariant(chromosome="chr1", position=1510035, reference="GCC", alternative=["G"])
        self._assert_metrics(variant=variant, expected_ac={'G': 12}, expected_dp=13)

    def _assert_metrics(self, variant: VafatorVariant, expected_ac, expected_dp,
                        min_base_quality=0, min_mapping_quality=0):
        pileups = get_variant_pileup(
            variant=variant, bam=self.bam_reader,
            min_base_quality=min_base_quality, min_mapping_quality=min_mapping_quality)
        coverage_metrics = get_metrics(variant=variant, pileups=pileups)
        self.assertEqual(expected_ac, coverage_metrics.ac)
        self.assertEqual(expected_dp, coverage_metrics.dp)
