#!/bin/bash


source tests/assert.sh
test_data=`pwd`/vafator/tests/resources
output=output/workflow_test_02
mkdir -p $output

echo -e "sample_1\t"$test_data"/test_tumor_normal.vcf\t"$test_data"/TESTX_S1_L001.bam\t"$test_data"/TESTX_S1_L002.bam" > $test_data/test_input.txt
echo -e "sample_2\t"$test_data"/test_single_sample.vcf\t"$test_data"/TESTX_S1_L001.bam,"$test_data"/TESTX_S1_L002.bam\t"$test_data"/TESTX_S1_L001.bam,"$test_data"/TESTX_S1_L002.bam" >> $test_data/test_input.txt
nextflow main.nf -profile test --output $output --input_files $test_data/test_input.txt --skip_multiallelic_filter

test -s $output/sample_1/test_tumor_normal.vaf.vcf || { echo "Missing sample 1 output file!"; exit 1; }
test -s $output/sample_1/test_tumor_normal.vaf.filtered_multiallelics.vcf || { echo "Missing sample 1 output file!"; exit 1; }
test -s $output/sample_2/test_single_sample.vaf.vcf || { echo "Missing sample 2 output file!"; exit 1; }
test -s $output/sample_2/test_single_sample.vaf.filtered_multiallelics.vcf || { echo "Missing sample 1 output file!"; exit 1; }
