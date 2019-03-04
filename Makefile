all: build clean prepare start

build:
	docker build -t scilab .

clean:
	docker stop $(shell docker ps -aq) || true && docker rm $(shell docker ps -aq) || true
	rm -rf Results
	rm -rf outputs
	rm -f output_python

prepare:
	mkdir Results
	mkdir outputs

start:
	python3 execDockers.py &> output_python

archive:
	tar -cf Results-`date +'%Y%m%d'`.tar ./Results

extract:
	tar -xf *.tar
