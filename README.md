# VAFator

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5565744.svg)](https://doi.org/10.5281/zenodo.5565744)
[![PyPI version](https://badge.fury.io/py/vafator.svg)](https://badge.fury.io/py/vafator)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/vafator/badges/version.svg)](https://anaconda.org/bioconda/vafator)
[![Run unit tests](https://github.com/TRON-Bioinformatics/vafator/actions/workflows/unit_tests.yml/badge.svg?branch=master)](https://github.com/TRON-Bioinformatics/vafator/actions/workflows/unit_tests.yml)
[![License](https://img.shields.io/badge/license-MIT-green)](https://opensource.org/licenses/MIT)

Annotate the variants in a VCF file with technical annotations from one or more BAM files.

Install from PyPI (`pip install vafator`) or from bioconda (`conda install bioconda::vafator`). 

Annotations:

* **Allele frequency (AF)**: ratio of reads supporting the alternate allele. When more than one alternate allele is present then one value per alternate allele is provided.
* **Allele count (AC)**: count of reads supporting the alternate allele.  When more than one alternate allele is present then one value per alternate allele is provided.
* **Depth of coverage (DP)**: number of reads covering the position of the variant

Outputs a VCF with the following annotations in the INFO field for tumor and normal:
```
chr1    12345       .       A       G       .       PASS  tumor_af=0.0;tumor_ac=0;tumor_dp=89;normal_af=0.0196;normal_ac=1;normal_dp=51
chr2    12345       .       A       G,T       .       PASS  tumor_af=0.2,0.2;tumor_ac=2,2;tumor_dp=10;normal_af=0.0,0.0;normal_ac=0,0;normal_dp=10
```

Both tumor and normal BAMs are optional, it can annotate only with the tumor BAM and viceversa.

If more than one BAM is provided for either the tumor or the normal then the annotations are calculated across all BAMs 
and for also each of them separately (eg: `tumor_af` provides the allele frequency across all tumor BAMs, `tumor_af_1` 
and `tumor_af_2` provide the allele frequency on the first and second BAM respectively).


Run it as follows:
```
$ vafator --help
usage: vafator [-h] --input-vcf INPUT_VCF --output-vcf OUTPUT_VCF
               [--normal-bams NORMAL_BAMS [NORMAL_BAMS ...]]
               [--tumor-bams TUMOR_BAMS [TUMOR_BAMS ...]]
               [--mapping-quality MAPPING_QUALITY]
               [--base-call-quality BASE_CALL_QUALITY]
               [--prefix BASE_CALL_QUALITY]

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
  --prefix PREFIX
                        When provided the annotations are preceded by this prefix, otherwise the annotations
                        are named as tumor_af, normal_af, tumor_ac, normal_ac, tumor_dp and normal_dp

Copyright (c) 2019-2021 TRON gGmbH (See LICENSE for licensing details) 
```
