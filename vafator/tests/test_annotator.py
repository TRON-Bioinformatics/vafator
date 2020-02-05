import os
import pkg_resources
from unittest import TestCase
from cyvcf2 import VCF
from vafator.annotator import Annotator
import vafator.tests.utils as test_utils


class TestAnnotator(TestCase):

    def test_annotator(self):
        input_file = pkg_resources.resource_filename(__name__, "resources/test1.vcf")
        output_vcf = pkg_resources.resource_filename(__name__, "resources/results/test_annotator1_output.vcf")
        bam1 = pkg_resources.resource_filename(__name__, "resources/COLO_829_n1.bam")
        bam2 = pkg_resources.resource_filename(__name__, "resources/COLO_829_t1.bam")
        annotator = Annotator(input_vcf=input_file, output_vcf=output_vcf, normal_bams=[bam1], tumor_bams=[bam2])
        annotator.run()

        self.assertTrue(os.path.exists(output_vcf))
        n_variants_input = test_utils._get_count_variants(input_file)
        n_variants_output = test_utils._get_count_variants(output_vcf)
        self.assertTrue(n_variants_input == n_variants_output)

    def _get_info_at(self, input_file, chromosome, position, annotation):
        vcf = VCF(input_file)
        self.assertIsNotNone(vcf)
        for v in vcf:
            if v.POS == position and v.CHROM == chromosome:
                vcf.close()
                return v.INFO.get(annotation)
        vcf.close()
        return {}