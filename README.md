# VAFator

[![PyPI version](https://badge.fury.io/py/vafator.svg)](https://badge.fury.io/py/vafator)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/vafator/badges/version.svg)](https://anaconda.org/bioconda/vafator)
[![Run unit tests](https://github.com/TRON-Bioinformatics/vafator/actions/workflows/unit_tests.yml/badge.svg?branch=master)](https://github.com/TRON-Bioinformatics/vafator/actions/workflows/unit_tests.yml)
[![License](https://img.shields.io/badge/license-MIT-green)](https://opensource.org/licenses/MIT)

Annotate the variants in a VCF file with allele frequencies, allele counts and depth of coverage for alternate alleles 
from one or more BAM files.

Outputs a VCF with the following annotations in the INFO field for tumor and normal:
```
chr1    12345       .       A       G       .       PASS  tumor_af=0.0;tumor_ac=0;tumor_dp=89;normal_af=0.0196078431372549;normal_ac=1;normal_dp=51
```

Both tumor and normal BAMs are optional, it can annotate only with the tumor BAM and viceversa.
If more than one BAM is provided for either the tumor or the normal the frequencies and counts are added up.

Run it as follows:
```
$ vafator --help
usage: vafator [-h] --input-vcf INPUT_VCF --output-vcf OUTPUT_VCF
               [--normal-bams NORMAL_BAMS [NORMAL_BAMS ...]]
               [--tumor-bams TUMOR_BAMS [TUMOR_BAMS ...]]
               [--mapping-quality MAPPING_QUALITY]
               [--base-call-quality BASE_CALL_QUALITY]

vafator v0.1.0

optional arguments:
  -h, --help            show this help message and exit
  --input-vcf INPUT_VCF
                        The VCF to annotate (default: None)
  --output-vcf OUTPUT_VCF
                        The annotated VCF (default: None)
  --normal-bams NORMAL_BAMS [NORMAL_BAMS ...]
                        Whitespace-separated list of normal BAMs to analyse
                        (default: [])
  --tumor-bams TUMOR_BAMS [TUMOR_BAMS ...]
                        Whitespace-separated list of tumor BAMs to analyse
                        (default: [])
  --mapping-quality MAPPING_QUALITY
                        All reads with a mapping quality lower or equal than
                        this threshold will be filtered out (default: 0)
  --base-call-quality BASE_CALL_QUALITY
                        All bases with a base call quality lower or equal than
                        this threshold will be filtered out (default: 29)

Copyright (c) 2015-2020 TRON gGmbH (See LICENSE for licensing details) 
```
