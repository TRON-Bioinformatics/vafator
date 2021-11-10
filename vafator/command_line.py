#!/usr/bin/env python
import argparse
import sys
import logging
import vafator
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
    parser.add_argument("--normal-bams", dest="normal_bams", nargs='+', action="store", default=[],
                        help="Whitespace-separated list of normal BAMs to analyse")
    parser.add_argument("--tumor-bams", dest="tumor_bams", nargs='+', action="store", default=[],
                        help="Whitespace-separated list of tumor BAMs to analyse")
    parser.add_argument("--mapping-quality", dest="mapping_quality", action="store", type=int, default=1,
                        help="All reads with a mapping quality below this threshold will be filtered out")
    parser.add_argument("--base-call-quality", dest="base_call_quality", action="store", type=int, default=30,
                        help="All bases with a base call quality below this threshold will be filtered out")
    parser.add_argument("--prefix", dest="prefix", action="store", type=str, default=None,
                        help="When provided the annotations are preceded by this prefix, otherwise the annotations"
                             "are named as tumor_af, normal_af, tumor_ac, normal_ac, tumor_dp and normal_dp")

    args = parser.parse_args()

    logging.info("Vafator starting...")
    try:
        annotator = Annotator(
            input_vcf=args.input_vcf,
            output_vcf=args.output_vcf,
            normal_bams=args.normal_bams,
            tumor_bams=args.tumor_bams,
            mapping_qual_thr=args.mapping_quality,
            base_call_qual_thr=args.base_call_quality,
            prefix=args.prefix
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
    parser.add_argument("-N", "--max_CN", required=False, default=6, type=int,
                        help="maximum total copy number for each observed clone")
    args = parser.parse_args()
    run_vafator2decifer(args)
