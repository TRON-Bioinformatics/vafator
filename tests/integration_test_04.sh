#!/bin/bash


source tests/assert.sh
test_data=vafator/tests/resources
output=output/integration_test_04
input_vcf=$test_data/project.NIST.hc.snps.indels.chr1_1000000_2000000.vcf
output_vcf=$output/vafator.vcf

mkdir -p $output
vafator --input-vcf $input_vcf \
--output-vcf $output_vcf \
--bam my_normal $test_data/project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.chr1_1000000_2000000.bam \
--bam my_tumor $test_data/project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.chr1_1000000_2000000.bam \
--bam my_metastasis $test_data/project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.chr1_1000000_2000000.bam \
--mapping-quality 40 --base-call-quality 40

vafator2decifer --vcf_file $output/vafator.vcf \
--samples my_normal,my_tumor \
--cna_file $test_data/best.seg.minimal.ucn \
--out_dir $output \
--min_depth 1 \
--min_alt_depth 1 \
--min_vaf 0.1


test -s $output/decifer.input.tsv || { echo "Missing decifer input output file!"; exit 1; }
test -s $output/decifer.purity.tsv || { echo "Missing decifer input output file!"; exit 1; }
test -s $output/vafator.vcf || { echo "Missing SNPs bed output file!"; exit 1; }
