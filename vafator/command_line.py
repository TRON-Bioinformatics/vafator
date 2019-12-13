#!/usr/bin/env python
import argparse
import sys
import logging
import vafator
from vafator.annotator import Annotator
from vafator.multiallelic_filter import MultiallelicFilter


epilog = "Copyright (c) 2015-2020 TRON gGmbH (See LICENSE for licensing details)"


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
    parser.add_argument("--mapping-quality", dest="mapping_quality", action="store", type=int, default=0,
                        help="All reads with a mapping quality lower or equal than this threshold will be filtered out")
    parser.add_argument("--base-call-quality", dest="base_call_quality", action="store", type=int, default=29,
                        help="All bases with a base call quality lower or equal than this threshold will be filtered out")

    args = parser.parse_args()

    logging.info("Vafator starting...")
    try:
        annotator = Annotator(
            input_vcf=args.input_vcf,
            output_vcf=args.output_vcf,
            normal_bams=args.normal_bams,
            tumor_bams=args.tumor_bams,
            mapping_qual_thr=args.mapping_quality,
            base_call_qual_thr=args.base_call_quality)
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
