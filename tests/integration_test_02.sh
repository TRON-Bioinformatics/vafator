#!/bin/bash


source tests/assert.sh
test_data=vafator/tests/resources
output=output/integration_test_02
input_vcf=$test_data/project.NIST.hc.snps.indels.chr1_1000000_2000000.vcf
output_vcf=$output/vafator.vcf
mkdir -p $output
vafator --input-vcf $input_vcf \
--output-vcf $output_vcf \
--normal-bams $test_data/project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.chr1_1000000_2000000.bam \
--tumor-bams $test_data/project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.chr1_1000000_2000000.bam

test -s $output/vafator.vcf || { echo "Missing VCF output file!"; exit 1; }
assert_eq `cat $output_vcf | grep -v '#' | grep normal_af | wc -l` `cat $input_vcf | grep -v '#' | wc -l` "Wrong number of variants"
assert_eq `cat $output_vcf | grep -v '#' | grep tumor_af | wc -l` `cat $input_vcf | grep -v '#' | wc -l` "Wrong number of variants"
assert_eq `cat $output_vcf | grep -v '#' | grep normal_ac | wc -l` `cat $input_vcf | grep -v '#' | wc -l` "Wrong number of variants"
assert_eq `cat $output_vcf | grep -v '#' | grep tumor_ac | wc -l` `cat $input_vcf | grep -v '#' | wc -l` "Wrong number of variants"
assert_eq `cat $output_vcf | grep -v '#' | grep normal_dp | wc -l` `cat $input_vcf | grep -v '#' | wc -l` "Wrong number of variants"
assert_eq `cat $output_vcf | grep -v '#' | grep tumor_dp | wc -l` `cat $input_vcf | grep -v '#' | wc -l` "Wrong number of variants"
