#!/usr/bin/env python
import argparse
import sys
import logging
import vafator
from vafator.power import DEFAULT_FPR, DEFAULT_ERROR_RATE
from vafator.hatchet2bed import run_hatchet2bed
from vafator.ploidies import PloidyManager
from vafator.annotator import Annotator
from vafator.multiallelic_filter import MultiallelicFilter
from vafator.vafator2decifer import run_vafator2decifer

epilog = "Copyright (c) 2019-2021 TRON gGmbH (See LICENSE for licensing details)"


def annotator():

    # set up logger
    parser = argparse.ArgumentParser(description="vafator v{}".format(vafator.VERSION),
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     epilog=epilog)
    parser.add_argument("--input-vcf", dest="input_vcf", action="store", help="The VCF to annotate", required=True)
    parser.add_argument("--output-vcf", dest="output_vcf", action="store", help="The annotated VCF", required=True)
    parser.add_argument('--bam', action='append', nargs=2,
                        metavar=('sample_name', 'bam_file'), default=[],
                        help='A sample name and a BAM file. Can be used multiple times to input multiple samples and '
                             'multiple BAM files. The same sample name can be used multiple times with different BAMs, '
                             'this will treated as replicates.')
    parser.add_argument("--mapping-quality", dest="mapping_quality", action="store", type=int, default=1,
                        help="All reads with a mapping quality below this threshold will be filtered out")
    parser.add_argument("--base-call-quality", dest="base_call_quality", action="store", type=int, default=30,
                        help="All bases with a base call quality below this threshold will be filtered out")
    parser.add_argument('--purity', action='append', nargs=2,
                        metavar=('sample_name', 'purity'), default=[],
                        help='A sample name and a tumor purity value. Can be used multiple times to input multiple '
                             'samples in combination with --bam. If no purity is provided for a given sample the '
                             'default value is 1.0')
    parser.add_argument("--tumor-ploidy", action='append', nargs=2,
                        metavar=('sample_name', 'tumor_ploidy'), default=[],
                        help='A sample name and a tumor ploidy. Can be used multiple times to input multiple '
                             'samples in combination with --bam. The tumor ploidy can be provided as a genome-wide '
                             'value (eg: --tumor-ploidy primary 2) or as local copy numbers in a BED file '
                             '(eg: --tumor-ploidy primary /path/to/copy_numbers.bed), see the documentation for '
                             'expected BED format (default: 2)')
    parser.add_argument("--normal-ploidy", dest="normal_ploidy", required=False, default=2, type=int,
                        help="Normal ploidy for the power calculation (default: 2)")
    parser.add_argument("--fpr", dest="fpr", required=False, default=DEFAULT_FPR, type=float,
                        help="False Positive Rate (FPR) to use in the power calculation")
    parser.add_argument("--error-rate", dest="error_rate", required=False, default=DEFAULT_ERROR_RATE, type=float,
                        help="Error rate to use in the power calculation")
    parser.add_argument("--include-ambiguous-bases", dest="include_ambiguous_bases", action='store_true',
                        help="Flag indicating to include ambiguous bases from the DP calculation")

    args = parser.parse_args()

    logging.info("Vafator starting...")

    bams = {}
    for sample_name, bam in args.bam:
        if sample_name in bams:
            bams[sample_name].append(bam)
        else:
            bams[sample_name] = [bam]

    purities = {}
    for sample_name, purity in args.purity:
        if sample_name in purities:
            raise ValueError('Multiple purity values provided for sample: {}'.format(sample_name))
        if sample_name not in bams:
            raise ValueError('Provided a purity value for a sample for which no BAM is provided: {}'.format(sample_name))
        purities[sample_name] = float(purity)

    tumor_ploidies = {}
    for sample_name, tumor_ploidy in args.tumor_ploidy:
        if sample_name in tumor_ploidies:
            raise ValueError('Multiple tumor ploidy values provided for sample: {}'.format(sample_name))
        if sample_name not in bams:
            raise ValueError(
                'Provided a tumor ploidy value for a sample for which no BAM is provided: {}'.format(sample_name))
        try:
            # checks if a genome-wide purity value was passed
            tumor_ploidies[sample_name] = PloidyManager(genome_wide_ploidy=float(tumor_ploidy))
        except ValueError:
            # checks if the non float-like value is a path to an existing file
            tumor_ploidies[sample_name] = PloidyManager(local_copy_numbers=tumor_ploidy)

    if len(bams) == 0:
        raise ValueError("Please, provide at least one bam file with '--bam sample_name /path/to/file.bam'")

    try:
        annotator = Annotator(
            input_vcf=args.input_vcf,
            output_vcf=args.output_vcf,
            input_bams=bams,
            mapping_qual_thr=args.mapping_quality,
            base_call_qual_thr=args.base_call_quality,
            purities=purities,
            tumor_ploidies=tumor_ploidies,
            normal_ploidy=int(args.normal_ploidy),
            fpr=args.fpr,
            error_rate=args.error_rate,
            include_ambiguous_bases=args.include_ambiguous_bases
        )
        annotator.run()
    except Exception as e:
        logging.error(str(e))
        sys.exit(-1)
    logging.info("Vafator finished!")


def multiallelics_filter():

    # set up logger
    parser = argparse.ArgumentParser(description="vafator v{}".format(vafator.VERSION),
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter, epilog=epilog)
    parser.add_argument("--input-vcf", dest="input_vcf", action="store", help="The VCF to annotate", required=True)
    parser.add_argument("--output-vcf", dest="output_vcf", action="store", help="The annotated VCF", required=True)
    parser.add_argument("--tumor-sample-name", dest="tumor_sample_name", action="store",
                        help='The tumor sample name (will look for annotation ${SAMPLE_NAME}_af)', default='tumor')
    args = parser.parse_args()

    logging.info("Vafator multiallelic filter starting...")
    try:
        filter = MultiallelicFilter(
            input_vcf=args.input_vcf,
            output_vcf=args.output_vcf,
            tumor_sample_name=args.tumor_sample_name
        )
        filter.run()
    except Exception as e:
        logging.error(str(e))
        sys.exit(-1)
    logging.info("Vafator multiallelic filter finished!")


def vafator2decifer():
    parser = argparse.ArgumentParser(description='Generate input for Decifer using VCF file and HATCHet CNA file')
    parser.add_argument("-V", "--vcf_file", required=True, type=str, help="single or multi-sample VCF file")
    parser.add_argument("-S", "--samples", required=True, type=str,
                        help="comma separated list of sample name prefixes to use for VAFator annotations, "
                             "eg: primary_tumor,metastasis_tumor; the annotations primary_tumor_ac, primary_tumor_dp, "
                             "etc. will be expected to exist")
    parser.add_argument("-C", "--cna_file", required=True, type=str, help="HATCHet CNA file: best.seg.ucn ")
    parser.add_argument("-O", "--out_dir", required=True, default="./", type=str,
                        help="directory for printing files; please make unique for each patient!")
    parser.add_argument("-M", "--min_depth", required=True, type=int, help="minimum depth PER sample")
    parser.add_argument("-A", "--min_alt_depth", required=True, type=int,
                        help="minimum depth of ALT allele in at least one sample")
    parser.add_argument("-F", "--min_vaf", required=True, type=float,
                        help="minimum VAF of ALT allele in at least one sample")
    parser.add_argument("-N", "--max_CN", required=False, default=6, type=int,
                        help="maximum total copy number for each observed clone")
    parser.add_argument("-B", "--exclude_list", required=False, default=None, type=str,
                        help="BED file of genomic regions to exclude")
    parser.add_argument("-p", "--min_purity", required=False, default=0.0, type=float,
                        help="minimum purity to consider samples")
    parser.add_argument("--snp_file", required=False, default=None, type=str,
                        help="HATCHet file containing germline SNP counts in tumor samples, baf/tumor.1bed")
    args = parser.parse_args()
    run_vafator2decifer(args)


def hatchet2bed():
    parser = argparse.ArgumentParser(description='Generate input for Decifer using VCF file and HATCHet CNA file')
    parser.add_argument("-i", "--input-file", required=True, type=str, help="input *.ucn hatchet file")
    parser.add_argument("-o", "--output-prefix", required=True, type=str,
                        help="output BED file prefix, one file will be created per sample in the input with the "
                             "average tumor copy number in each segment")
    args = parser.parse_args()
    run_hatchet2bed(input_file=args.input_file, output_prefix=args.output_prefix)
