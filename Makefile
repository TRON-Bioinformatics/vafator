
all : clean test check

clean:
	rm -rf output
	#rm -rf work
	rm -f .nextflow.log*
	rm -rf .nextflow*


test:
	nextflow main.nf --help
	echo "sample_1\t"`pwd`"/test_data/test_tumor_normal.vcf\t"`pwd`"/test_data/TESTX_S1_L001.bam\t"`pwd`"/test_data/TESTX_S1_L002.bam\n" > test_data/test_input.txt
	echo "sample_2\t"`pwd`"/test_data/test_single_sample.vcf\t"`pwd`"/test_data/TESTX_S1_L001.bam,"`pwd`"/test_data/TESTX_S1_L002.bam\t"`pwd`"/test_data/TESTX_S1_L001.bam,"`pwd`"/test_data/TESTX_S1_L002.bam" >> test_data/test_input.txt
	nextflow main.nf -profile test,conda --output output/test1 --input_files test_data/test_input.txt
	nextflow main.nf -profile test,conda --output output/test2 --input_files test_data/test_input.txt --skip_multiallelic_filter
	nextflow main.nf -profile test,conda --output output/test3 --input_files test_data/test_input.txt --prefix RNA --skip_multiallelic_filter


check:
	test -s output/test1/sample_1/test_tumor_normal.vaf.vcf || { echo "Missing test 1 sample 1 output file!"; exit 1; }
	test -s output/test1/sample_1/test_tumor_normal.vaf.filtered_multiallelics.vcf || { echo "Missing test 1 sample 1 output file!"; exit 1; }
	test -s output/test1/sample_2/test_single_sample.vaf.vcf || { echo "Missing test 1 sample 2 output file!"; exit 1; }
	test -s output/test1/sample_2/test_single_sample.vaf.filtered_multiallelics.vcf || { echo "Missing test 1 sample 1 output file!"; exit 1; }
	test -s output/test2/sample_1/test_tumor_normal.vaf.vcf || { echo "Missing test 2 sample 1 output file!"; exit 1; }
	test -s output/test2/sample_2/test_single_sample.vaf.vcf || { echo "Missing test 2 sample 2 output file!"; exit 1; }
	test -s output/test3/sample_1/test_tumor_normal.vaf.vcf || { echo "Missing test 3 sample 1 output file!"; exit 1; }
	test -s output/test3/sample_2/test_single_sample.vaf.vcf || { echo "Missing test 3 sample 2 output file!"; exit 1; }

integration_tests:
	bash tests/integration_test_00.sh
	bash tests/integration_test_01.sh
	bash tests/integration_test_02.sh
	bash tests/integration_test_03.sh


