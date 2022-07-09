from unittest import TestCase

import pkg_resources

from vafator.hatchet2bed import run_hatchet2bed

class Hatchet2bedTest(TestCase):

    def hatchet2bed_test(self):
        run_hatchet2bed(
            input_file=pkg_resources.resource_filename(__name__, "resources/best.seg.minimal.ucn"))