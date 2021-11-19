
all : clean workflow_tests

clean:
	rm -rf output
	#rm -rf work
	rm -f .nextflow.log*
	rm -rf .nextflow*


workflow_tests:
	bash tests/workflow_test_00.sh
	bash tests/workflow_test_01.sh
	bash tests/workflow_test_02.sh

integration_tests:
	bash tests/integration_test_00.sh
	bash tests/integration_test_01.sh
	bash tests/integration_test_02.sh
	bash tests/integration_test_03.sh


