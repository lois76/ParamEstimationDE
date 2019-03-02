all: build clean prepare start

build:
	docker build -t scilab .

clean:
	rm -rf Results
	rm -rf outputs
	rm output_python

prepare:
	mkdir Results
	mkdir outputs

start:
	python3 execDockers.py > output_python

archive:
	tar -cvf Results-`date +'%Y%m%d'`.tar ./Results

extract:
	tar -xvf *.tar