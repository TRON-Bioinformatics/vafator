from collections import Counter
from unittest import TestCase
from unittest.mock import MagicMock
import pkg_resources
import pysam

from vafator.tests.utils import VafatorVariant
from vafator.pileups import (
    get_variant_pileup, _get_metrics_from_column,
    _get_insertion_metrics_from_column, _get_deletion_metrics_from_column,
)
from vafator.pileups import VariantRecord


def _make_pileup_col(reads):
    col = MagicMock()
    col.pileups = reads
    return col


class FakeAlignment:
    def __init__(self, reference_start, cigartuples, mapping_quality, query=None,
                 query_sequence=None, query_qualities=None):
        self.reference_start = reference_start
        self.cigartuples = cigartuples
        self.mapping_quality = mapping_quality
        self.query = query
        self.query_sequence = query_sequence
        self.query_qualities = query_qualities


def _make_insertion_read(indel, reference_start, cigartuples, query, mapping_quality=60,
                         query_position_or_next=0):
    """Build a mock pileup read with an insertion (indel > 0)."""
    read = MagicMock()
    read.indel = indel
    read.alignment = FakeAlignment(
        reference_start=reference_start,
        cigartuples=cigartuples,
        query=query,
        mapping_quality=mapping_quality,
    )
    read.query_position_or_next = query_position_or_next
    return read


def _make_ref_read(mapping_quality=60, query_position_or_next=0):
    """Build a mock pileup read with no indel (indel == 0), counted as reference."""
    read = MagicMock()
    read.indel = 0
    read.alignment = FakeAlignment(
        reference_start=0, cigartuples=[], mapping_quality=mapping_quality)
    read.query_position_or_next = query_position_or_next
    return read


def _make_deletion_read(indel, reference_start, cigartuples, mapping_quality=60,
                        query_position_or_next=0):
    """Build a mock pileup read with a deletion (indel < 0)."""
    read = MagicMock()
    read.indel = indel
    read.alignment = FakeAlignment(
        reference_start=reference_start,
        cigartuples=cigartuples,
        mapping_quality=mapping_quality,
    )
    read.query_position_or_next = query_position_or_next
    return read


class TestInsertionMetricsFromColumn(TestCase):
    # CIGAR op codes: 0=M, 1=I, 2=D, 3=N, 4=S, 7==, 8=X

    def test_no_reads(self):
        variant = VariantRecord(CHROM='chr1', POS=100, REF='A', ALT=['ATT'])
        col = _make_pileup_col([])
        metrics = _get_insertion_metrics_from_column(variant, col)
        self.assertEqual(metrics.dp, 0)
        self.assertEqual(metrics.ac['ATT'], 0)

    def test_one_matching_insertion(self):
        # variant: A -> ATT at pos 100 (insertion of "TT", length 2)
        variant = VariantRecord(CHROM='chr1', POS=100, REF='A', ALT=['ATT'])
        read = _make_insertion_read(
            indel=2,
            reference_start=99,
            cigartuples=[(0, 1), (1, 2)],
            query='ATT',
            mapping_quality=60,
            query_position_or_next=1
        )
        col = _make_pileup_col([read])
        metrics = _get_insertion_metrics_from_column(variant, col)
        self.assertEqual(metrics.dp, 1)
        self.assertEqual(metrics.ac['ATT'], 1)
        self.assertEqual(metrics.mqs['ATT'], 60.0)

    def test_insertion_wrong_sequence_not_counted(self):
        variant = VariantRecord(CHROM='chr1', POS=100, REF='A', ALT=['ATT'])
        read = _make_insertion_read(
            indel=2,
            reference_start=99,
            cigartuples=[(0, 1), (1, 2)],
            query='ATC',
            mapping_quality=60,
            query_position_or_next=1
        )
        col = _make_pileup_col([read])
        metrics = _get_insertion_metrics_from_column(variant, col)
        self.assertEqual(metrics.dp, 1)
        self.assertEqual(metrics.ac['ATT'], 0)

    def test_insertion_wrong_length_not_counted(self):
        variant = VariantRecord(CHROM='chr1', POS=100, REF='A', ALT=['ATT'])
        read = _make_insertion_read(
            indel=3,
            reference_start=99,
            cigartuples=[(0, 1), (1, 3)],
            query='ATTT',
            mapping_quality=60,
            query_position_or_next=1
        )
        col = _make_pileup_col([read])
        metrics = _get_insertion_metrics_from_column(variant, col)
        self.assertEqual(metrics.ac['ATT'], 0)

    def test_ref_reads_counted_in_mq(self):
        variant = VariantRecord(CHROM='chr1', POS=100, REF='A', ALT=['ATT'])
        ref_read = _make_ref_read(mapping_quality=55, query_position_or_next=30)
        col = _make_pileup_col([ref_read])
        metrics = _get_insertion_metrics_from_column(variant, col)
        self.assertEqual(metrics.dp, 1)
        self.assertEqual(metrics.ac['ATT'], 0)
        self.assertEqual(metrics.mqs['A'], 55.0)

    def test_mixed_insertion_and_ref_reads(self):
        variant = VariantRecord(CHROM='chr1', POS=100, REF='A', ALT=['ATT'])
        ins_read = _make_insertion_read(
            indel=2, reference_start=99,
            cigartuples=[(0, 1), (1, 2)],
            query='ATT', mapping_quality=60, query_position_or_next=1
        )
        ref_read = _make_ref_read(mapping_quality=55)
        col = _make_pileup_col([ins_read, ref_read])
        metrics = _get_insertion_metrics_from_column(variant, col)
        self.assertEqual(metrics.dp, 2)
        self.assertEqual(metrics.ac['ATT'], 1)


class TestDeletionMetricsFromColumn(TestCase):
    # CIGAR op codes: 0=M, 2=D, 3=N, 7==, 8=X

    def test_no_reads(self):
        variant = VariantRecord(CHROM='chr1', POS=100, REF='ATT', ALT=['A'])
        col = _make_pileup_col([])
        metrics = _get_deletion_metrics_from_column(variant, col)
        self.assertEqual(metrics.dp, 0)
        self.assertEqual(metrics.ac['A'], 0)

    def test_one_matching_deletion(self):
        # variant: ATT -> A at pos 100 (deletion of length 2)
        # read starts at 98, M2 advances to 100, D2 matches the deletion
        variant = VariantRecord(CHROM='chr1', POS=100, REF='ATT', ALT=['A'])
        read = _make_deletion_read(
            indel=-2,
            reference_start=98,
            cigartuples=[(0, 2), (2, 2)],
            mapping_quality=60,
            query_position_or_next=2
        )
        col = _make_pileup_col([read])
        metrics = _get_deletion_metrics_from_column(variant, col)
        self.assertEqual(metrics.dp, 1)
        self.assertEqual(metrics.ac['A'], 1)
        self.assertEqual(metrics.mqs['A'], 60.0)

    def test_deletion_wrong_length_not_counted(self):
        # deletion of length 3 but variant expects length 2
        variant = VariantRecord(CHROM='chr1', POS=100, REF='ATT', ALT=['A'])
        read = _make_deletion_read(
            indel=-3,
            reference_start=98,
            cigartuples=[(0, 2), (2, 3)],
            mapping_quality=60,
            query_position_or_next=2
        )
        col = _make_pileup_col([read])
        metrics = _get_deletion_metrics_from_column(variant, col)
        self.assertEqual(metrics.ac['A'], 0)

    def test_deletion_wrong_position_not_counted(self):
        # deletion starts at 99 not 100
        variant = VariantRecord(CHROM='chr1', POS=100, REF='ATT', ALT=['A'])
        read = _make_deletion_read(
            indel=-2,
            reference_start=98,
            cigartuples=[(0, 1), (2, 2)],  # M1 advances to 99, not 100
            mapping_quality=60,
            query_position_or_next=1
        )
        col = _make_pileup_col([read])
        metrics = _get_deletion_metrics_from_column(variant, col)
        self.assertEqual(metrics.ac['A'], 0)

    def test_ref_reads_counted_in_mq(self):
        variant = VariantRecord(CHROM='chr1', POS=100, REF='ATT', ALT=['A'])
        ref_read = _make_ref_read(mapping_quality=45, query_position_or_next=10)
        col = _make_pileup_col([ref_read])
        metrics = _get_deletion_metrics_from_column(variant, col)
        self.assertEqual(metrics.dp, 1)
        self.assertEqual(metrics.ac['A'], 0)
        self.assertEqual(metrics.mqs['ATT'], 45.0)

    def test_mixed_deletion_and_ref_reads(self):
        variant = VariantRecord(CHROM='chr1', POS=100, REF='ATT', ALT=['A'])
        del_read = _make_deletion_read(
            indel=-2, reference_start=98,
            cigartuples=[(0, 2), (2, 2)],
            mapping_quality=60, query_position_or_next=2
        )
        ref_read = _make_ref_read(mapping_quality=55)
        col = _make_pileup_col([del_read, ref_read])
        metrics = _get_deletion_metrics_from_column(variant, col)
        self.assertEqual(metrics.dp, 2)
        self.assertEqual(metrics.ac['A'], 1)


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

    def _assert_metrics(self, variant, expected_ac, expected_dp,
                        min_base_quality=0, min_mapping_quality=0):
        pileups = get_variant_pileup(
            variant=variant, bam=self.bam_reader,
            min_base_quality=min_base_quality, min_mapping_quality=min_mapping_quality)
        pileup_col = next(iter(pileups), None)
        if pileup_col is None:
            coverage_metrics = None
        else:
            coverage_metrics = _get_metrics_from_column(
                variant=variant, pileup_col=pileup_col, include_ambiguous_bases=True)
        if expected_dp == 0:
            self.assertIsNone(coverage_metrics)
        else:
            self.assertEqual(expected_ac, coverage_metrics.ac)
            self.assertEqual(expected_dp, coverage_metrics.dp)