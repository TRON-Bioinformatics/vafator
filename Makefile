
all : clean integration_tests

clean:
	rm -rf output
	#rm -rf work
	rm -f .nextflow.log*
	rm -rf .nextflow*

integration_tests:
	bash tests/integration_test_00.sh
	bash tests/integration_test_01.sh
	bash tests/integration_test_02.sh
	bash tests/integration_test_03.sh
