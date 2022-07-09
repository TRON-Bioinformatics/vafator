<img src="images/VAFator-logos.jpeg" width="200">

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5565743.svg)](https://doi.org/10.5281/zenodo.5565743)
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
* **Power**: probability that a given somatic mutation is undetected given the total coverage, supporting reads and expected allele frequency.

VAFator uses cyvcf2 (Pederson, 2017) to read/write VCF files and pysam (https://github.com/pysam-developers/pysam) to read BAM files.
Both libraries are cython wrappers around HTSlib (Bonfield, 2021).

## How to install

Install from PyPI (`pip install vafator`) or from bioconda (`conda install bioconda::vafator`).

When installaing from PyPI there are some system dependencies that will need to be met: 
* libcurl
* libz
* liblzma
* htslib=1.14
* 

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
`normal_dp`, `normal_af`, `normal_pw`, `primary_ac`, `primary_dp`, etc. 

## Support for technical replicates

If more than one BAM  for the same sample is provided then the annotations are calculated across all BAMs 
and for also each of them separately (eg: `primary_af` provides the allele frequency across all primary tumor BAMs, 
`primary_af_1` and `primary_af_2` provide the allele frequency on the first and second BAM respectively).

```
vafator --input-vcf /path/yo/your.vcf \
--output-vcf /path/yo/your_vafator.vcf \ 
--bam primary /path/to/your_primary_tumor_1.bam \
--bam primary /path/to/your_primary_tumor_2.bam
```

### Read filtering

Use the parameters `--mapping-quality` and `--base-call-quality` to define the minimum quality values for each read.
All reads with quality values below these thresholds will be filtered out. 
The default values are MQ >= 1 and BQ >= 30.

Overlapping reads from read pairs are not double counted. The read with the highest base call quality is chosen.

Reads flagged as duplicates are not counted.


### Power (or probability of an undetected somatic mutation)

We estimate the probability that there is an undetected somatic mutation at a specific genomic location, given the
- n: observed number of reads at the position (coverage)
- k: observed number of reads supporting the mutation (variant reads)
- f: expected variant allele frequency (VAF)

We model this with a binomial distribution `binom(n, f, k)`.

The expected VAF is by default 0.5, making several assumptions: 
1) no normal contamination in tumor sample (default purity: 1.0)
2) ploidy in normal of 2
3) ploidy in tumor of 2
4) the mutation is clonal and the mutation multiplicity is one irrespective of the copy number

Expected VAF (f) is calculated as follows:

![Expected VAF](images/expected_vaf_formula.png)

The purity in a given tumor sample can be specified with the parameter `--purity primary 0.5`.

The ploidy in the normal sample can be specified as `--normal_ploidy 3`. 
Only one normal sample may be specified.
Also, in the normal only genome-wide ploidy may be specified.

The ploidy in a given tumor sample can be specified either as a genome-wide value `--tumor-ploidy metastasis 4`
or local copy numbers in a BED file as follows `--tumor-ploidy metastasis /path/to/copy_numbers.bed`.

The expected format of the BedGraph file ([specification](https://genome.ucsc.edu/goldenpath/help/bedgraph.html)) 
without any track information as indicated below:
```
chr1    10000   20000   3.2
chr1    20000   30000   2.6
[...]
```

**NOTE**: beware that unlike VCF files where genomic positions are 1-based, BedGraph positions are 0-based. The intervals
in a BED file are half-closed.
**NOTE 2**: local copy numbers may be float numbers
**NOTE 3**: beware that no multiple test correction is applied to the power

In order to integrate Hatchet copy numbers were multiple clones are considered, there is a specific 
command `hatchet2bed` that transforms a Hatchet .ucn output file into one BedGraph file as described above per sample.
The Hatchet reported clone copy numbers are averaged by the relative abundance of each clone. 

Run this as follows:
```
hatchet2bed --input-file your.hatchet.ucn --output-prefix your_beds_prefix
```

The Hatchet file expects the following columns, where any number of clones is supported:
```
#CHR	START	END	SAMPLE	cn_normal	u_normal	cn_clone1	u_clone1	cn_clone2	u_clone2	cn_clone3	u_clone3	cn_clone4	u_clone4
```

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
nextflow run tron-bioinformatics/tronflow-vcf-postprocessing -r 2.2.0 -profile conda ...
```

See https://github.com/TRON-Bioinformatics/tronflow-vcf-postprocessing for more details


## Support for indels

VAFator provides equivalent annotations for indels. Depth of coverage and allele frequency are calculated on the 
position immediately before the indel. Only insertions and deletions as recorded in the CIGAR matching the respective 
coordinates and sequence from the VCF file are taken into account. Any read supporting a similar but not identical indel
is not counted. 

**NOTE**: multiallelic mutations are not supported for indels, the indel in the multiallelic position will be 
annotated with null values.

## Support for MNVs

Not supported at the moment when not decomposed.

If running the nextflow pipeline indicated above, MNVs and complex variants are by default decomposed and hence
correctly annotated by VAFator.

## Bibliography

- Pedersen, B. S., & Quinlan, A. R. (2017). cyvcf2: fast, flexible variant analysis with Python. Bioinformatics, 33(12), 1867–1869. https://doi.org/10.1093/BIOINFORMATICS/BTX057
- Bonfield, J. K., Marshall, J., Danecek, P., Li, H., Ohan, V., Whitwham, A., Keane, T., & Davies, R. M. (2021). HTSlib: C library for reading/writing high-throughput sequencing data. GigaScience, 10(2). https://doi.org/10.1093/GIGASCIENCE/GIAB007
