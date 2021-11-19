# VAFator

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5565744.svg)](https://doi.org/10.5281/zenodo.5565744)
[![PyPI version](https://badge.fury.io/py/vafator.svg)](https://badge.fury.io/py/vafator)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/vafator/badges/version.svg)](https://anaconda.org/bioconda/vafator)
[![Run unit tests](https://github.com/TRON-Bioinformatics/vafator/actions/workflows/unit_tests.yml/badge.svg?branch=master)](https://github.com/TRON-Bioinformatics/vafator/actions/workflows/unit_tests.yml)
[![License](https://img.shields.io/badge/license-MIT-green)](https://opensource.org/licenses/MIT)

VAFator annotates the variants in a VCF file with technical annotations from multiple BAM files. 
Supports annotating somatic variant calls with the annotations from the normal and the tumor samples; although
it can also be used for germline variant calls.

Annotations:

* **Allele frequency (AF)**: ratio of reads supporting the alternate allele.
* **Allele count (AC)**: count of reads supporting the alternate allele. 
* **Depth of coverage (DP)**: number of reads covering the position of the variant


## How to install

Install from PyPI (`pip install vafator`) or from bioconda (`conda install bioconda::vafator`).

When installaing from PyPI there are some system dependencies that will need to be met: 
* libcurl
* libz
* liblzma
* htslib=1.14

## How to run

Run as follows:
```
vafator --input-vcf /path/yo/your.vcf \
--output-vcf /path/yo/your_vafator.vcf \ 
--bam normal /path/to/your_normal.bam \
--bam primary /path/to/your_primary_tumor.bam \
--bam metastasis /path/to/your_metastasis_tumor.bam
```

This will add annotations for each of the three samples `normal`, `primary` and `metastasis`: `normal_ac`, 
`normal_dp`, `normal_af`, `primary_ac`, `primary_dp`, `primary_af`, 
`metastasis_ac`, `metastasis_dp` and `metastasis_af`. 

If more than one BAM  for the same sample is provided then the annotations are calculated across all BAMs 
and for also each of them separately (eg: `primary_af` provides the allele frequency across all primary tumor BAMs, 
`primary_af_1` and `primary_af_2` provide the allele frequency on the first and second BAM respectively).

```
vafator --input-vcf /path/yo/your.vcf \
--output-vcf /path/yo/your_vafator.vcf \ 
--bam primary /path/to/your_primary_tumor_1.bam \
--bam primary /path/to/your_primary_tumor_2.bam
```

Alternatively, you can use `--normal-bams` and/or `--tumor-bams` and the sample names will be predefined to `normal` 
and `tumor`respectively.

```
vafator --input-vcf /path/yo/your.vcf \
--output-vcf /path/yo/your_vafator.vcf \ 
--normal-bams /path/to/your_normal.bam \
--tumor-bams /path/to/your_tumor_1.bam,/path/to/your_tumor_2.bam
```

Use the parameters `--mapping-quality` and `--base-call-quality` to define the minimum quality values for each read.
All reads with quality values velow these thresholds will be filtered out.

Overlapping reads from read pairs are not double counted. The read with the highest base call quality is chosen.

## Understanding the output

The output is a VCF with the some new annotations in the INFO field for the provided sample names.
The example below contains vafator annotations for two samples named `normal` and `tumor`.
```
chr1    12345       .       A       G       .       PASS  tumor_af=0.0;tumor_ac=0;tumor_dp=89;normal_af=0.0196;normal_ac=1;normal_dp=51
chr2    12345       .       A       G,T       .       PASS  tumor_af=0.2,0.2;tumor_ac=2,2;tumor_dp=10;normal_af=0.0,0.0;normal_ac=0,0;normal_dp=10
```

**NOTE**: notice that VAFator does not annotate samples in the FORMAT field, but in the INFO field

## Filter for multi-allelic variants

Multi-allelic variants are those that have more than one alternative allele (e.g.: A>C,G).
This tool allows to select the allele with the highest allele frequency and filter out the lower frequency allele.

Run as follows:
```
multiallelics-filter --input-vcf /path/to/your_vafator.vcf \
--output-vcf /path/to/your_vafator_filtered.vcf \
--tumor-sample-name <SAMPLE> 
```

The above will look for the annotation `<SAMPLE>_af` and for multi-allelic variants it will filter out those with lower 
frequencies. Beware, that if the multi-allelic variants are split into more than one line in the VCF nothing will be 
filtered out.

## Run as a Nextflow pipeline

VAFator is available as a Nextflow pipeline for convenience.

Run as follows:
```
nextflow run tron-bioinformatics/vafator -r 1.1.0 -profile conda --input_files /path/to/your.tsv
```

where `--input_files` expects three tab-separated columns **without a header**:

| Patient name     | VCF                    | Sample name:BAM              |
|-----------------|------------------------|-------------------------|
| patient_1        | /path/to/sample_1.vcf  | my_tumor:/path/to/sample_1_tumor_1.bam  |
| patient_1        | /path/to/sample_1.vcf  | my_tumor:/path/to/sample_1_tumor_2.bam   |
| patient_1        | /path/to/sample_1.vcf  | my_metastasis:/path/to/sample_1_metastasis_1.bam  |
| patient_1        | /path/to/sample_1.vcf  | my_normal:/path/to/sample_1_normal.bam   |
| patient_2        | /path/to/sample_2.vcf  | my_tumor:/path/to/sample_2_tumor.bam   |
| patient_2        | /path/to/sample_2.vcf  | my_normal:/path/to/sample_1_normal.bam   |

The above corresponds to two patients, `patient_1` and `patient_2`. Where for `patient_1` we have four samples of which
two replicates for a sample named `my_tumor` and two additional samples named `my_metastasis` and `my_normal`.

Optional parameters:

* `--output`: the folder where to publish output
* `--skip_multiallelic_filter`: skip the filtering of multiallelics by frequency in the tumor (only highest frequency 
  variant at the same position is kept)
* `--base_call_quality`: minimum base call quality, reads with values below will be filtered out (default: 30)
* `--mapping_quality`: minimum mapping quality, reads with values below will be filtered out (default: 1)


## Support for indels

VAFator provides equivalent annotations for indels. Depth of coverage and allele frequency are calculated on the 
position immediately before the indel. Only insertions and deletions as recorded in the CIGAR matching the respective 
coordinates and sequence from the VCF file are taken into account. Any read supporting a similar but not equal indel
will not be counted. 
Also, multiallelic mutations are not supported for indels.


## Support for MNVs

Not supported
