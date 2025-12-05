.PHONY: build run clean

build:
	docker build -t msa-benchmark .

run:
	docker run -v $(PWD)/data:/app/data -v $(PWD)/results:/app/results msa-benchmark

clean:
	docker system prune -f
	rm -rf results/*

all: build run
