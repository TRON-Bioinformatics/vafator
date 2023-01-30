import os

import numpy as np
import pkg_resources
from unittest import TestCase
from cyvcf2 import VCF

from vafator.ploidies import PloidyManager
from vafator.annotator import Annotator
import vafator.tests.utils as test_utils
import time
from logzero import logger


EXPECTED_ANNOTATIONS = [
    'af', 'dp', 'ac', 'n', 'pu', 'pw', 'k', 'eaf', 'bq', 'mq', 'pos', 'rsmq', 'rsmq_pv', 'rsbq', 'rsbq_pv', 'rspos',
    'rspos_pv'
]

# replicates do not have EAF annotation
EXPECTED_ANNOTATIONS_REPLICATES = [
    'af', 'dp', 'ac', 'n', 'pu', 'pw', 'k', 'bq', 'mq', 'pos', 'rsmq', 'rsmq_pv', 'rsbq', 'rsbq_pv', 'rspos',
    'rspos_pv'
]


class TestAnnotator(TestCase):

    def test_annotator(self):
        input_file = pkg_resources.resource_filename(__name__, "resources/test1.vcf")
        output_vcf = pkg_resources.resource_filename(__name__, "resources/results/test_annotator1_output.vcf")
        bam1 = pkg_resources.resource_filename(__name__, "resources/COLO_829_n1.bam")
        bam2 = pkg_resources.resource_filename(__name__, "resources/COLO_829_t1.bam")
        annotator = Annotator(
            input_vcf=input_file, output_vcf=output_vcf, input_bams={"normal": [bam1], "tumor": [bam2]})
        annotator.run()

        self.assertTrue(os.path.exists(output_vcf))
        n_variants_input = test_utils._get_count_variants(input_file)
        n_variants_output = test_utils._get_count_variants(output_vcf)
        self.assertTrue(n_variants_input == n_variants_output)

        info_annotations = test_utils._get_info_fields(output_vcf)
        for a in EXPECTED_ANNOTATIONS:
            self.assertTrue("tumor_{}".format(a) in info_annotations)
            self.assertTrue("normal_{}".format(a) in info_annotations)

    def test_annotator_with_multiple_bams(self):
        input_file = pkg_resources.resource_filename(__name__, "resources/test1.vcf")
        output_vcf = pkg_resources.resource_filename(__name__, "resources/results/test_annotator1_output.vcf")
        bam1 = pkg_resources.resource_filename(__name__, "resources/COLO_829_n1.bam")
        bam2 = pkg_resources.resource_filename(__name__, "resources/COLO_829_t1.bam")
        annotator = Annotator(
            input_vcf=input_file, output_vcf=output_vcf, input_bams={"normal": [bam1, bam2], "tumor": [bam1, bam2]})
        annotator.run()

        self.assertTrue(os.path.exists(output_vcf))
        n_variants_input = test_utils._get_count_variants(input_file)
        n_variants_output = test_utils._get_count_variants(output_vcf)
        self.assertTrue(n_variants_input == n_variants_output)

        info_annotations = test_utils._get_info_fields(output_vcf)
        for a in EXPECTED_ANNOTATIONS_REPLICATES:
            self.assertTrue("tumor_{}_1".format(a) in info_annotations,
                            "Missing annotation tumor_{}_1".format(a))
            self.assertTrue("normal_{}_1".format(a) in info_annotations,
                            "Missing annotation normal_{}_1".format(a))

    def test_annotator_with_prefix(self):
        input_file = pkg_resources.resource_filename(__name__, "resources/test1.vcf")
        output_vcf = pkg_resources.resource_filename(__name__, "resources/results/test_annotator1_output.vcf")
        bam1 = pkg_resources.resource_filename(__name__, "resources/COLO_829_n1.bam")
        bam2 = pkg_resources.resource_filename(__name__, "resources/COLO_829_t1.bam")
        annotator = Annotator(
            input_vcf=input_file, output_vcf=output_vcf,
            input_bams={"RNA_normal": [bam1, bam2], "RNA_tumor": [bam1, bam2]})
        annotator.run()

        self.assertTrue(os.path.exists(output_vcf))
        n_variants_input = test_utils._get_count_variants(input_file)
        n_variants_output = test_utils._get_count_variants(output_vcf)
        self.assertTrue(n_variants_input == n_variants_output)

        info_annotations = test_utils._get_info_fields(output_vcf)
        for a in EXPECTED_ANNOTATIONS_REPLICATES:
            self.assertTrue("RNA_tumor_{}_1".format(a) in info_annotations,
                            "Missing annotation RNA_tumor_{}_1".format(a))
            self.assertTrue("RNA_normal_{}_1".format(a) in info_annotations,
                            "Missing annotation RNA_normal_{}_1".format(a))

    def test_annotator_with_mnvs(self):
        input_file = pkg_resources.resource_filename(__name__, "resources/test_tumor_normal.vcf")
        output_vcf = pkg_resources.resource_filename(__name__, "resources/results/test_tumor_normal_output.vcf")
        bam1 = pkg_resources.resource_filename(__name__, "resources/COLO_829_n1.bam")
        bam2 = pkg_resources.resource_filename(__name__, "resources/COLO_829_t1.bam")
        annotator = Annotator(
            input_vcf=input_file, output_vcf=output_vcf,
            input_bams={"RNA_normal": [bam1, bam2], "RNA_tumor": [bam1, bam2]})
        annotator.run()

        self.assertTrue(os.path.exists(output_vcf))
        n_variants_input = test_utils._get_count_variants(input_file)
        n_variants_output = test_utils._get_count_variants(output_vcf)
        self.assertTrue(n_variants_input == n_variants_output)

        info_annotations = test_utils._get_info_fields(output_vcf)
        for a in EXPECTED_ANNOTATIONS_REPLICATES:
            self.assertTrue("RNA_tumor_{}_1".format(a) in info_annotations,
                            "Missing annotation RNA_tumor_{}_1".format(a))
            self.assertTrue("RNA_normal_{}_1".format(a) in info_annotations,
                            "Missing annotation RNA_normal_{}_1".format(a))

    def _get_info_at(self, input_file, chromosome, position, annotation):
        vcf = VCF(input_file)
        self.assertIsNotNone(vcf)
        for v in vcf:
            if v.POS == position and v.CHROM == chromosome:
                vcf.close()
                return v.INFO.get(annotation)
        vcf.close()
        return {}

    def test_nist(self):
        input_file = pkg_resources.resource_filename(
            __name__, "resources/project.NIST.hc.snps.indels.chr1_1000000_2000000.vcf")
        output_vcf = pkg_resources.resource_filename(
            __name__, "resources/results/project.NIST.hc.snps.indels.chr1_1000000_2000000.vaf.vcf")
        bam_file = pkg_resources.resource_filename(
            __name__,
            "resources/project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.chr1_1000000_2000000.bam")
        start = time.time()
        annotator = Annotator(input_vcf=input_file, output_vcf=output_vcf, input_bams={"normal": [bam_file]})
        annotator.run()
        duration = time.time() - start
        logger.info("Duration {} seconds".format(round(duration, 3)))

        self._assert_vafator_vcf(output_vcf, sample_name='normal')

        n_variants_input = test_utils._get_count_variants(input_file)
        n_variants_output = test_utils._get_count_variants(output_vcf)
        self.assertTrue(n_variants_input == n_variants_output)

        info_annotations = test_utils._get_info_fields(output_vcf)
        self.assertTrue("normal_af" in info_annotations)
        self.assertTrue("normal_ac" in info_annotations)
        self.assertTrue("normal_dp" in info_annotations)

        variant = test_utils._get_mutation_at_position(output_vcf, 'chr1', 1506035)
        self.assertIsNotNone(variant)
        self.assertEqual(variant.INFO['normal_ac'], 3)
        self.assertEqual(variant.INFO['normal_dp'], 3)
        self.assertEqual(variant.INFO['normal_af'], 1.0)
        self.assertEqual(variant.INFO['normal_pu'], 1.0)
        self.assertEqual(variant.INFO['normal_eaf'], 0.5)
        self.assertEqual(variant.INFO['normal_mq'][0], 0)
        self.assertEqual(variant.INFO['normal_mq'][1], 60.0)
        self.assertEqual(variant.INFO['normal_bq'][0], 0)
        self.assertEqual(variant.INFO['normal_bq'][1], 37.0)
        self.assertEqual(variant.INFO['normal_pos'][0], 0)
        self.assertEqual(variant.INFO['normal_pos'][1], 56.0)

        variant = test_utils._get_mutation_at_position(output_vcf, 'chr1', 1509825)
        self.assertIsNotNone(variant)
        self.assertEqual(variant.INFO['normal_ac'], 13)
        self.assertEqual(variant.INFO['normal_dp'], 20)
        # these values are rounded to six digits inside the VCF, not sure why when read the representation is
        # different...
        self.assertEqual(round(variant.INFO['normal_af'], 5), 0.65)
        self.assertEqual(round(variant.INFO['normal_pu'], 5), 0.94234)
        self.assertEqual(variant.INFO['normal_eaf'], 0.5)
        self.assertEqual(variant.INFO['normal_mq'][0], 60.0)
        self.assertEqual(variant.INFO['normal_mq'][1], 60.0)
        self.assertEqual(variant.INFO['normal_bq'][0], 37.0)
        self.assertEqual(variant.INFO['normal_bq'][1], 35.0)
        self.assertEqual(variant.INFO['normal_pos'][0], 31.0)
        self.assertEqual(variant.INFO['normal_pos'][1], 41.0)

        # this is a deletion
        variant = test_utils._get_mutation_at_position(output_vcf, 'chr1', 1323143)
        self.assertIsNotNone(variant)
        self.assertEqual(variant.INFO['normal_ac'], 20)
        self.assertEqual(variant.INFO['normal_dp'], 21)
        # these values are rounded to six digits inside the VCF, not sure why when read the representation is
        # different...
        self.assertEqual(round(variant.INFO['normal_af'], 5), 0.95238)
        self.assertEqual(round(variant.INFO['normal_pu'], 5), 1.0)
        self.assertEqual(variant.INFO['normal_eaf'], 0.5)
        self.assertEqual(variant.INFO['normal_mq'][0], 60.0)
        self.assertEqual(variant.INFO['normal_mq'][1], 60.0)
        self.assertEqual(variant.INFO['normal_bq'][0], 0.0)
        self.assertEqual(variant.INFO['normal_bq'][1], 0.0)
        self.assertEqual(variant.INFO['normal_pos'][0], 50.0)
        self.assertEqual(variant.INFO['normal_pos'][1], 21.0)

        # this is an insertion
        variant = test_utils._get_mutation_at_position(output_vcf, 'chr1', 1935367)
        self.assertIsNotNone(variant)
        self.assertEqual(variant.INFO['normal_ac'], 1)
        self.assertEqual(variant.INFO['normal_dp'], 2)
        # these values are rounded to six digits inside the VCF, not sure why when read the representation is
        # different...
        self.assertEqual(round(variant.INFO['normal_af'], 5), 0.5)
        self.assertEqual(round(variant.INFO['normal_pu'], 5), 0.75)
        self.assertEqual(variant.INFO['normal_eaf'], 0.5)
        self.assertEqual(variant.INFO['normal_mq'][0], 60.0)
        self.assertEqual(variant.INFO['normal_mq'][1], 29.0)
        self.assertEqual(variant.INFO['normal_bq'][0], 0.0)
        self.assertEqual(variant.INFO['normal_bq'][1], 0.0)
        self.assertEqual(variant.INFO['normal_pos'][0], 95.0)
        self.assertEqual(variant.INFO['normal_pos'][1], 56.0)

    def test_nist_with_replicates(self):
        input_file = pkg_resources.resource_filename(
            __name__, "resources/project.NIST.hc.snps.indels.chr1_1000000_2000000.vcf")
        output_vcf = pkg_resources.resource_filename(
            __name__, "resources/results/project.NIST.hc.snps.indels.chr1_1000000_2000000.vaf_replicates.vcf")
        bam_file = pkg_resources.resource_filename(
            __name__,
            "resources/project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.chr1_1000000_2000000.bam")
        start = time.time()
        annotator = Annotator(input_vcf=input_file, output_vcf=output_vcf, input_bams={"normal": [bam_file]})
        annotator.run()
        duration = time.time() - start
        logger.info("Duration {} seconds".format(round(duration, 3)))

        self._assert_vafator_vcf(output_vcf, sample_name='normal')
        self._assert_vafator_vcf(output_vcf, sample_name='normal', replicate=1)
        self._assert_vafator_vcf(output_vcf, sample_name='normal', replicate=2)

    def test_annotator_bams_order(self):
        input_file = pkg_resources.resource_filename(__name__, "resources/test1.vcf")
        output_vcf = pkg_resources.resource_filename(__name__, "resources/results/test_annotator1_output.vcf")
        output_vcf_2 = pkg_resources.resource_filename(__name__, "resources/results/test_annotator2_output.vcf")
        bam1 = pkg_resources.resource_filename(__name__, "resources/COLO_829_n1.bam")
        bam2 = pkg_resources.resource_filename(__name__, "resources/COLO_829_t1.bam")

        Annotator(input_vcf=input_file, output_vcf=output_vcf, input_bams={"normal": [bam1], "tumor": [bam2]}).run()
        Annotator(input_vcf=input_file, output_vcf=output_vcf_2, input_bams={"tumor": [bam2], "normal": [bam1]}).run()

        self.assertTrue(os.path.exists(output_vcf))
        self.assertTrue(os.path.exists(output_vcf_2))

        self._assert_vafator_vcf(output_vcf, sample_name='normal')
        self._assert_vafator_vcf(output_vcf, sample_name='tumor')
        self._assert_vafator_vcf(output_vcf_2, sample_name='normal')
        self._assert_vafator_vcf(output_vcf_2, sample_name='tumor')

        vcf = VCF(output_vcf)
        vcf_2 = VCF(output_vcf_2)

        for v, v2 in zip(vcf, vcf_2):
            for a in EXPECTED_ANNOTATIONS:
                self.assertEqual(
                    v.INFO.get("normal_{}".format(a), ""),
                    v2.INFO.get("normal_{}".format(a), ""),
                    "Variant {}:{}:{}>{} is missing annotation normal_{}".format(v.CHROM, v.POS, v.REF, v.ALT[0], a))
                self.assertEqual(
                    v.INFO.get("tumor_{}".format(a), ""),
                    v2.INFO.get("tumor_{}".format(a), ""),
                    "Variant {}:{}:{}>{} is missing annotation tumor_{}".format(v.CHROM, v.POS, v.REF, v.ALT[0], a))

    def test_annotator_with_purities(self):
        input_file = pkg_resources.resource_filename(__name__, "resources/test1.vcf")
        output_vcf = pkg_resources.resource_filename(__name__, "resources/results/test_annotator1_output.vcf")
        bam1 = pkg_resources.resource_filename(__name__, "resources/COLO_829_n1.bam")
        bam2 = pkg_resources.resource_filename(__name__, "resources/COLO_829_t1.bam")
        annotator = Annotator(
            input_vcf=input_file, output_vcf=output_vcf, input_bams={"normal": [bam1], "tumor": [bam2]},
            purities={"tumor": 0.8}, tumor_ploidies={"tumor": PloidyManager(genome_wide_ploidy=2.8)}
        )
        annotator.run()

        self.assertTrue(os.path.exists(output_vcf))
        n_variants_input = test_utils._get_count_variants(input_file)
        n_variants_output = test_utils._get_count_variants(output_vcf)
        self.assertTrue(n_variants_input == n_variants_output)

        info_annotations = test_utils._get_info_fields(output_vcf)
        for a in EXPECTED_ANNOTATIONS:
            self.assertTrue("tumor_{}".format(a) in info_annotations)
            self.assertTrue("normal_{}".format(a) in info_annotations)

        annotator = Annotator(
            input_vcf=input_file, output_vcf=output_vcf, input_bams={"normal": [bam1], "tumor": [bam2]},
            purities={"tumor": 0.2}, tumor_ploidies={"tumor": PloidyManager(genome_wide_ploidy=1.5)}
        )
        annotator.run()

    def _assert_vafator_vcf(self, vcf_filename, sample_name, replicate=None):
        self.assertTrue(os.path.exists(vcf_filename))
        vcf = VCF(vcf_filename)
        for v in vcf:
            # p-values or VAFs
            self._assert_probability(v.INFO.get(self._get_annotation_name('rsbq_pv', sample_name, replicate=replicate), 0))
            self._assert_probability(v.INFO.get(self._get_annotation_name('rsmq_pv', sample_name, replicate=replicate), 0))
            self._assert_probability(v.INFO.get(self._get_annotation_name('rspos_pv', sample_name, replicate=replicate), 0))
            self._assert_probability(v.INFO.get(self._get_annotation_name('eaf', sample_name, replicate=replicate), 0))
            self._assert_probability(v.INFO.get(self._get_annotation_name('af', sample_name, replicate=replicate), 0))

            # positive integer annotations
            self._assert_positive_integer(v.INFO.get(self._get_annotation_name('ac', sample_name, replicate=replicate), 0))
            self._assert_positive_integer(v.INFO.get(self._get_annotation_name('dp', sample_name, replicate=replicate), 0))
            self._assert_positive_integer(v.INFO.get(self._get_annotation_name('mq', sample_name, replicate=replicate), 0))
            self._assert_positive_integer(v.INFO.get(self._get_annotation_name('bq', sample_name, replicate=replicate), 0))
            self._assert_positive_integer(v.INFO.get(self._get_annotation_name('pos', sample_name, replicate=replicate), 0))
        vcf.close()

    @staticmethod
    def _get_annotation_name(annotation_name, sample_name, replicate=None):
        annotation_name = "{}_{}".format(sample_name, annotation_name)
        if replicate:
            annotation_name = "{}_{}".format(annotation_name, replicate)
        return annotation_name

    def _assert_probability(self, annotation):
        if isinstance(annotation, list) or isinstance(annotation, tuple):
            for a in annotation:
                self.assertTrue(0.0 <= float(a) <= 1.0, "Expected probability has a value of {}".format(a))
        else:
            self.assertTrue(0.0 <= float(annotation) <= 1.0, "Expected probability has a value of {}".format(annotation))

    def _assert_positive_integer(self, annotation):
        if isinstance(annotation, list) or isinstance(annotation, tuple):
            for a in annotation:
                self.assertTrue(np.isnan(a) or 0.0 <= a, "Expected positive integer has a value of {}".format(a))
        else:
            self.assertTrue(np.isnan(annotation) or 0.0 <= annotation, "Expected positive integer has a value of {}".format(annotation))
