.PHONY: build run clean

# Build the Docker image
build:
	docker build -t msa-benchmark .

# Run the benchmarking
run:
	docker run -v $(PWD)/data:/app/data -v $(PWD)/results:/app/results msa-benchmark

# Clean up
clean:
	docker system prune -f
	rm -rf results/*

# Run everything in sequence
all: build run