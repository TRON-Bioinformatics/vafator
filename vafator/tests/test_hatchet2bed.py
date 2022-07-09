import os
from unittest import TestCase
import pkg_resources
from vafator.hatchet2bed import run_hatchet2bed


class Hatchet2bedTest(TestCase):

    def test_hatchet2bed(self):
        run_hatchet2bed(
            input_file=pkg_resources.resource_filename(__name__, "resources/best.seg.minimal.ucn"),
            output_prefix=pkg_resources.resource_filename(__name__, "resources/best.seg.minimal")
        )
        self.assertTrue(
            os.path.exists(pkg_resources.resource_filename(__name__, "resources/best.seg.minimal.my_tumor.bed")))
        self.assertTrue(
            os.path.exists(pkg_resources.resource_filename(__name__, "resources/best.seg.minimal.my_metastasis.bed")))
