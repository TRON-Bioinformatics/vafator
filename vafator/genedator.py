import gtfparse
import pysam
import numpy as np
import csv
import argparse
from logzero import logger
import sys


DEFAULT_BIOTYPES = ["protein_coding", "IG_V_gene", "IG_C_gene", "IG_J_gene", "IG_D_gene",
                    "TR_C_gene", "TR_J_gene", "TR_V_gene", "TR_D_gene", ]
BATCH_SIZE = 10000


class Genedator:

    def __init__(self, gff_file, bam_file, output_csv_file,  biotypes=DEFAULT_BIOTYPES):
        gtf = gtfparse.read_gtf(gff_file)
        self.genes = gtf[(gtf.feature == "gene") & (gtf.gene_type.isin(biotypes))]
        self.bam = pysam.AlignmentFile(bam_file, "rb")
        self.fd = open(output_csv_file, 'w', newline='')
        self.writer = csv.writer(self.fd, delimiter='\t')

    def run(self):
        data = []
        counter = 0
        for i, g in self.genes.iterrows():
            chromosome = g.seqname
            gene_name = g.gene_name
            for read in self.bam.fetch(chromosome, g.start, g.end):
                if not read.is_duplicate:
                    data.append((chromosome, gene_name, read.query_name,
                                 read.alen,
                                 read.cigarstring,
                                 read.flag,
                                 read.isize,
                                 read.mapping_quality,
                                 np.median(read.query_qualities),
                                 np.max(read.query_qualities),
                                 np.min(read.query_qualities)))
                    if len(data) > BATCH_SIZE:
                        self.writer.writerows(data)
                        counter = counter + len(data)
                        data = []
                        logger.info("Wrote {} reads".format(counter))
        self.fd.close()


if __name__ == '__main__':
    # set up logger
    parser = argparse.ArgumentParser(
        description="genedator", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input-bam", dest="input_bam", action="store", help="The BAM file to annotate with", required=True)
    parser.add_argument("--output-csv", dest="output_csv", action="store", help="The annotations CSV", required=True)
    parser.add_argument("--gff", dest="gff", action="store",
                        help='The GFF file defining the genes to annotate', default='tumor')
    args = parser.parse_args()

    try:
        logger.info("Starting Genedator!")
        Genedator(
            bam_file=args.input_bam,
            output_csv_file=args.output_csv,
            gff_file=args.gff
        ).run()
    except Exception as e:
        logger.exception(e)
        sys.exit(-1)
    logger.info("Genedator finished!")
