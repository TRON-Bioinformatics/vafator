#!/bin/bash


source tests/assert.sh
test_data=`pwd`/vafator/tests/resources
output=output/workflow_test_01
mkdir -p $output

echo -e "sample_1\t"$test_data"/test_tumor_normal.vcf\tmy_tumor:"$test_data"/TESTX_S1_L001.bam" > $test_data/test_input.txt
echo -e "sample_1\t"$test_data"/test_tumor_normal.vcf\tmy_normal:"$test_data"/TESTX_S1_L002.bam" >> $test_data/test_input.txt
echo -e "sample_2\t"$test_data"/test_single_sample.vcf\tmy_tumor:"$test_data"/TESTX_S1_L001.bam" >> $test_data/test_input.txt
echo -e "sample_2\t"$test_data"/test_single_sample.vcf\tmy_tumor:"$test_data"/TESTX_S1_L002.bam" >> $test_data/test_input.txt
echo -e "sample_2\t"$test_data"/test_single_sample.vcf\tmy_normal:"$test_data"/TESTX_S1_L001.bam" >> $test_data/test_input.txt
nextflow main.nf -profile test --output $output --input_files $test_data/test_input.txt

test -s $output/sample_1/test_tumor_normal.vaf.vcf || { echo "Missing sample 1 output file!"; exit 1; }
test -s $output/sample_1/test_tumor_normal.vaf.filtered_multiallelics.vcf || { echo "Missing sample 1 output file!"; exit 1; }
test -s $output/sample_2/test_single_sample.vaf.vcf || { echo "Missing sample 2 output file!"; exit 1; }
test -s $output/sample_2/test_single_sample.vaf.filtered_multiallelics.vcf || { echo "Missing sample 1 output file!"; exit 1; }

assert_eq `cat $output/sample_1/test_tumor_normal.vaf.vcf | grep -v '#' | grep my_normal_af | wc -l` `cat $test_data/test_tumor_normal.vcf | grep -v '#' | wc -l` "Wrong number of variants"
assert_eq `cat $output/sample_1/test_tumor_normal.vaf.vcf | grep -v '#' | grep my_tumor_af | wc -l` `cat $test_data/test_tumor_normal.vcf | grep -v '#' | wc -l` "Wrong number of variants"

