from unittest import TestCase

import pkg_resources
import pysam
from logzero import logger


class TestAnnotator(TestCase):

    def test_insertions_pileups(self):
        bam_file = pkg_resources.resource_filename(
            __name__, "resources/project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.chr1_0_1000000.bam")
        bam_reader = pysam.AlignmentFile(bam_file)
        chromosome = "chr1"
        position = 866511
        length = 4

        ac = 0
        dp = 0

        pileups = bam_reader.pileup(contig=chromosome, start=position - 1, stop=position, truncate=True)
        pileup = next(pileups)
        logger.info("Pileup at position {}".format(pileup.pos))
        logger.info("Number of reads {}".format(len(pileup.pileups)))
        if pileup.pileups:
            logger.info("indel={}".format(sum([r.indel for r in pileup.pileups])))
            logger.info("is_del={}".format(sum([r.is_del for r in pileup.pileups])))
            for r in pileup.pileups:
                dp += 1
                if r.indel > 0 and r.is_del == 0:
                    # read with an insertion
                    start = r.alignment.reference_start
                    for c, l in r.alignment.cigartuples:
                        if c in [0, 2]:
                            start += l
                            if start > position:
                                break
                        elif c == 1:
                            if start == position and l == length:
                                ac +=1
            logger.info("AC={}".format(ac))
            logger.info("DP={}".format(dp))
            af = float(ac) / dp
            logger.info("AF={}".format(af))
            self.assertEqual(ac, 15)
            self.assertEqual(dp, 25)

    def test_deletions_pileups_drafting(self):
        bam_file = pkg_resources.resource_filename(
            __name__, "resources/project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.chr1_1000000_2000000.bam")
        bam_reader = pysam.AlignmentFile(bam_file)
        chromosome = "chr1"
        pileups = bam_reader.pileup(contig=chromosome, start=1887090, stop=1887091, truncate=True)
        is_del = 0
        indel = 0
        for pileup in pileups:
            if pileup.pileups:
                indel += sum([r.indel for r in pileup.pileups])
                is_del += sum([r.is_del for r in pileup.pileups])
                pass
        logger.info("indel={}".format(indel))
        logger.info("is_del={}".format(is_del))

    def test_deletions_pileups(self):
        bam_file = pkg_resources.resource_filename(
            __name__, "resources/project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.chr1_0_1000000.bam")
        bam_reader = pysam.AlignmentFile(bam_file)
        chromosome = "chr1"
        position = 970561
        length = 2

        ac = 0
        dp = 0

        pileups = bam_reader.pileup(contig=chromosome, start=position - 1, stop=position, truncate=True)
        pileup = next(pileups)
        logger.info("Pileup at position {}".format(pileup.pos))
        logger.info("Number of reads {}".format(len(pileup.pileups)))
        if pileup.pileups:
            logger.info("indel={}".format(sum([r.indel for r in pileup.pileups])))
            logger.info("is_del={}".format(sum([r.is_del for r in pileup.pileups])))
            for r in pileup.pileups:
                dp += 1
                if r.indel > 0 and r.is_del == 0:
                    # read with an insertion
                    start = r.alignment.reference_start
                    for c, l in r.alignment.cigartuples:
                        if c in [0]:
                            start += l
                            if start > position:
                                break
                        elif c == 2:
                            if start == position and l == length:
                                ac +=1
                            else:
                                start += l
                        if start > position:
                            break
            logger.info("AC={}".format(ac))
            logger.info("DP={}".format(dp))
            af = float(ac) / dp
            logger.info("AF={}".format(af))
            self.assertEqual(ac, 15)
            self.assertEqual(dp, 25)


