
all : clean integration_tests

clean:
	rm -rf output
	#rm -rf work
	rm -f .nextflow.log*
	rm -rf .nextflow*

integration_tests:
	bash tests/integration_test_00.sh
	bash tests/integration_test_01.sh
	bash tests/integration_test_03.sh
	bash tests/integration_test_04.sh
	bash tests/integration_test_05.sh
	bash tests/integration_test_06.sh

jupyter:
	export USER_ID=$$(id -u $$USER); export USER_GROUP=$$(id -g $$USER); docker-compose up -d --build

jupyter-down:
	docker-compose down
