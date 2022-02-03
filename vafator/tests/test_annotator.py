import os
import pkg_resources
from unittest import TestCase
from cyvcf2 import VCF
from vafator.annotator import Annotator
import vafator.tests.utils as test_utils
import time
from logzero import logger


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
        self.assertTrue("tumor_af" in info_annotations)
        self.assertTrue("normal_af" in info_annotations)
        self.assertTrue("tumor_ac" in info_annotations)
        self.assertTrue("normal_ac" in info_annotations)
        self.assertTrue("tumor_dp" in info_annotations)
        self.assertTrue("normal_dp" in info_annotations)

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
        self.assertTrue("tumor_af_1" in info_annotations)
        self.assertTrue("normal_af_1" in info_annotations)
        self.assertTrue("tumor_ac_1" in info_annotations)
        self.assertTrue("normal_ac_1" in info_annotations)
        self.assertTrue("tumor_dp_1" in info_annotations)
        self.assertTrue("normal_dp_1" in info_annotations)
        self.assertTrue("tumor_af_2" in info_annotations)
        self.assertTrue("normal_af_2" in info_annotations)
        self.assertTrue("tumor_ac_2" in info_annotations)
        self.assertTrue("normal_ac_2" in info_annotations)
        self.assertTrue("tumor_dp_2" in info_annotations)
        self.assertTrue("normal_dp_2" in info_annotations)

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
        self.assertTrue("RNA_tumor_af_1" in info_annotations)
        self.assertTrue("RNA_normal_af_1" in info_annotations)
        self.assertTrue("RNA_tumor_ac_1" in info_annotations)
        self.assertTrue("RNA_normal_ac_1" in info_annotations)
        self.assertTrue("RNA_tumor_dp_1" in info_annotations)
        self.assertTrue("RNA_normal_dp_1" in info_annotations)
        self.assertTrue("RNA_tumor_af_2" in info_annotations)
        self.assertTrue("RNA_normal_af_2" in info_annotations)
        self.assertTrue("RNA_tumor_ac_2" in info_annotations)
        self.assertTrue("RNA_normal_ac_2" in info_annotations)
        self.assertTrue("RNA_tumor_dp_2" in info_annotations)
        self.assertTrue("RNA_normal_dp_2" in info_annotations)

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
        self.assertTrue("RNA_tumor_af_1" in info_annotations)
        self.assertTrue("RNA_normal_af_1" in info_annotations)
        self.assertTrue("RNA_tumor_ac_1" in info_annotations)
        self.assertTrue("RNA_normal_ac_1" in info_annotations)
        self.assertTrue("RNA_tumor_dp_1" in info_annotations)
        self.assertTrue("RNA_normal_dp_1" in info_annotations)
        self.assertTrue("RNA_tumor_af_2" in info_annotations)
        self.assertTrue("RNA_normal_af_2" in info_annotations)
        self.assertTrue("RNA_tumor_ac_2" in info_annotations)
        self.assertTrue("RNA_normal_ac_2" in info_annotations)
        self.assertTrue("RNA_tumor_dp_2" in info_annotations)
        self.assertTrue("RNA_normal_dp_2" in info_annotations)

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

        self.assertTrue(os.path.exists(output_vcf))
        n_variants_input = test_utils._get_count_variants(input_file)
        n_variants_output = test_utils._get_count_variants(output_vcf)
        self.assertTrue(n_variants_input == n_variants_output)

        info_annotations = test_utils._get_info_fields(output_vcf)
        self.assertTrue("normal_af" in info_annotations)
        self.assertTrue("normal_ac" in info_annotations)
        self.assertTrue("normal_dp" in info_annotations)

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

        vcf = VCF(output_vcf)
        vcf_2 = VCF(output_vcf_2)

        for v, v2 in zip(vcf, vcf_2):
            self.assertEqual(v.INFO["normal_dp"], v2.INFO["normal_dp"])
            self.assertEqual(v.INFO["tumor_dp"], v2.INFO["tumor_dp"])
