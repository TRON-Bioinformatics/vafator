import os
import pkg_resources
from unittest import TestCase
from cyvcf2 import VCF
from vafator.multiallelic_filter import MultiallelicFilter
import vafator.tests.utils as test_utils


class TestMultiallelicFilter(TestCase):

    def setUp(self):
        # self.reference_genome_reader = TestReferenceGenomeReader()
        pass

    def test_no_variants_filtered(self):
        input_file = pkg_resources.resource_filename(__name__, "resources/test1.vcf")
        output_vcf = pkg_resources.resource_filename(__name__, "resources/results/test1_output.vcf")
        multiallelic_filter = MultiallelicFilter(input_vcf=input_file, output_vcf=output_vcf, tumor_sample_name='tumor')
        multiallelic_filter.run()

        self.assertTrue(os.path.exists(output_vcf))
        n_variants_input = test_utils._get_count_variants(input_file)
        n_variants_output = test_utils._get_count_variants(output_vcf)
        self.assertTrue(n_variants_input == n_variants_output)

    def test_two_variants_filtered(self):
        input_file = pkg_resources.resource_filename(__name__, "resources/test2.vcf")
        output_vcf = pkg_resources.resource_filename(__name__, "resources/results/test2_output.vcf")
        multiallelic_filter = MultiallelicFilter(input_vcf=input_file, output_vcf=output_vcf, tumor_sample_name='tumor')
        multiallelic_filter.run()

        self.assertTrue(os.path.exists(output_vcf))
        n_variants_input = test_utils._get_count_variants(input_file)
        n_variants_output = test_utils._get_count_variants(output_vcf)
        self.assertTrue(n_variants_input == n_variants_output + 3, "input:{}; output:{}".format(
            n_variants_input, n_variants_output))

        af1 = self._get_info_at(output_vcf, chromosome="chr4", position=1235, annotation='tumor_af')
        self.assertTrue(af1, 0.2)
        multiallelic1 = self._get_info_at(output_vcf, chromosome="chr4", position=1235, annotation='multiallelic')
        self.assertTrue(multiallelic1, "T,0.1")

        af2 = self._get_info_at(output_vcf, chromosome="chr6", position=1235, annotation='tumor_af')
        self.assertTrue(af2, 0.2)
        multiallelic2 = self._get_info_at(output_vcf, chromosome="chr6", position=1235, annotation='multiallelic')
        self.assertTrue(multiallelic2, "A,0.01")

        af3 = self._get_info_at(output_vcf, chromosome="chr6", position=1234, annotation='tumor_af')
        self.assertTrue(af3, 0.5)

    def test_different_reference_is_kept(self):
        input_file = pkg_resources.resource_filename(__name__, "resources/test3.vcf")
        output_vcf = pkg_resources.resource_filename(__name__, "resources/results/test3_output.vcf")
        multiallelic_filter = MultiallelicFilter(input_vcf=input_file, output_vcf=output_vcf, tumor_sample_name='tumor')
        multiallelic_filter.run()

        self.assertTrue(os.path.exists(output_vcf))
        n_variants_input = test_utils._get_count_variants(input_file)
        n_variants_output = test_utils._get_count_variants(output_vcf)
        self.assertTrue(n_variants_input == n_variants_output)
        self.assertTrue(n_variants_input == 2)

    def test_three_multiallelics(self):
        input_file = pkg_resources.resource_filename(__name__, "resources/test4.vcf")
        output_vcf = pkg_resources.resource_filename(__name__, "resources/results/test4_output.vcf")
        multiallelic_filter = MultiallelicFilter(input_vcf=input_file, output_vcf=output_vcf, tumor_sample_name='tumor')
        multiallelic_filter.run()

        self.assertTrue(os.path.exists(output_vcf))
        n_variants_input = test_utils._get_count_variants(input_file)
        n_variants_output = test_utils._get_count_variants(output_vcf)
        self.assertTrue(n_variants_input == n_variants_output + 2)
        self.assertTrue(n_variants_output == 1)

    def test_equal_af(self):
        input_file = pkg_resources.resource_filename(__name__, "resources/test5.vcf")
        output_vcf = pkg_resources.resource_filename(__name__, "resources/results/test5_output.vcf")
        multiallelic_filter = MultiallelicFilter(input_vcf=input_file, output_vcf=output_vcf, tumor_sample_name='tumor')
        multiallelic_filter.run()

        self.assertTrue(os.path.exists(output_vcf))
        n_variants_input = test_utils._get_count_variants(input_file)
        n_variants_output = test_utils._get_count_variants(output_vcf)
        self.assertTrue(n_variants_input == n_variants_output + 2)
        self.assertTrue(n_variants_output == 1)

    def _get_info_at(self, input_file, chromosome, position, annotation):
        vcf = VCF(input_file)
        self.assertIsNotNone(vcf)
        for v in vcf:
            if v.POS == position and v.CHROM == chromosome:
                vcf.close()
                return v.INFO.get(annotation)
        vcf.close()
        return {}