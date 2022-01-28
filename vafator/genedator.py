import gtfparse
import pysam
import numpy as np
import csv
import argparse
from logzero import logger
import sys
import gzip


DEFAULT_BIOTYPES = ["protein_coding", "IG_V_gene", "IG_C_gene", "IG_J_gene", "IG_D_gene",
                    "TR_C_gene", "TR_J_gene", "TR_V_gene", "TR_D_gene", ]
BATCH_SIZE = 100000


class Genedator:

    def __init__(self, gtf_file, bam_file, output_csv_file, biotypes=DEFAULT_BIOTYPES):
        logger.info("Reading GTF file...")
        gtf = gtfparse.read_gtf(gtf_file)
        self.genes = gtf[(gtf.feature == "gene") & (gtf.gene_type.isin(biotypes))]
        self.bam = pysam.AlignmentFile(bam_file, "rb")
        self.fd = gzip.open(output_csv_file, 'wt')
        self.writer = csv.writer(self.fd, delimiter='\t')

    def run(self):
        # write header
        self.writer.writerow((
            "chromosome", "gene_name",
            "read_length_pct25", "read_length_pct50", "read_length_pct75",
            "mapping_quality_pct25", "mapping_quality_pct50", "mapping_quality_pct75",
            "insert_size_pct25", "insert_size_pct50", "insert_size_pct75",
            "median_bcq_pct25", "median_bcq_pct25", "median_bcq_pct25",
            "count_reads",
            "count_paired",
            "count_proper_pair",
            "count_qcfail",
            "count_read1",
            "count_read2",
            "count_reverse",
            "count_secondary",
            "count_supplementary",
            "count_duplicates"))
        num_genes = self.genes.shape[0]
        count_genes = 0
        for i, g in self.genes.iterrows():
            chromosome = g.seqname
            gene_name = g.gene_name
            read_lengths = []
            mqs = []
            insert_sizes = []
            median_bcq = []
            count_reads = 0
            count_duplicates = 0
            count_paired = 0
            count_proper_pair = 0
            count_qcfail = 0
            count_read1 = 0
            count_read2 = 0
            count_reverse = 0
            count_secondary = 0
            count_supplementary = 0
            for read in self.bam.fetch(chromosome, g.start, g.end):

                if not read.is_duplicate:
                    read_lengths.append(read.rlen)
                    mqs.append(read.mapping_quality)
                    insert_sizes.append(read.isize)
                    median_bcq.append(np.median(read.query_qualities))
                    count_reads += 1
                    if read.is_paired:
                        count_paired += 1
                    if read.is_proper_pair:
                        count_proper_pair += 1
                    if read.is_qcfail:
                        count_qcfail += 1
                    if read.is_read1:
                        count_read1 += 1
                    if read.is_read2:
                        count_read2 += 1
                    if read.is_reverse:
                        count_reverse += 1
                    if read.is_secondary:
                        count_secondary += 1
                    if read.is_supplementary:
                        count_supplementary += 1

                else:
                    count_duplicates += 1

            self.writer.writerow((
                chromosome, gene_name,
                np.percentile(read_lengths, 0.25),
                np.median(read_lengths),
                np.percentile(read_lengths, 0.75),
                np.percentile(mqs, 0.25),
                np.median(mqs),
                np.percentile(mqs, 0.75),
                np.percentile(insert_sizes, 0.25),
                np.median(insert_sizes),
                np.percentile(insert_sizes, 0.75),
                np.percentile(median_bcq, 0.25),
                np.median(median_bcq),
                np.percentile(median_bcq, 0.75),
                count_reads,
                count_paired,
                count_proper_pair,
                count_qcfail,
                count_read1,
                count_read2,
                count_reverse,
                count_secondary,
                count_supplementary,
                count_duplicates
              ))
            count_genes += 1
            if count_genes % 100 == 0:
                logger.info("Annotated {} % of genes".format(round(float(count_genes)/num_genes, 3)))
        self.fd.close()


if __name__ == '__main__':
    # set up logger
    parser = argparse.ArgumentParser(
        description="genedator", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input-bam", dest="input_bam", action="store", help="The BAM file to annotate with", required=True)
    parser.add_argument("--output-csv", dest="output_csv", action="store", help="The annotations CSV", required=True)
    parser.add_argument("--gtf", dest="gtf", action="store",
                        help='The GTF file defining the genes to annotate', default='tumor')
    args = parser.parse_args()

    try:
        logger.info("Starting Genedator!")
        Genedator(
            bam_file=args.input_bam,
            output_csv_file=args.output_csv,
            gtf_file=args.gtf
        ).run()
    except Exception as e:
        logger.exception(e)
        sys.exit(-1)
    logger.info("Genedator finished!")
