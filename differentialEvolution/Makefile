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

startParamEstimation:
	python3 execDockers.py &> output_python

startSeveralSameSimu:
	python3 execSeveralRuns.py &> output_python

archive:
	tar -czf Results-$(shell date +'%Y%m%d').tar.gz ./Results

extract:
	tar -zxf *.tar.gz
