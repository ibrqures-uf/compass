# COMPASS: MSA Benchmarking Tool

A benchmarking system for comparing Multiple Sequence Alignment (MSA) tools. We test five popular alignment programs (MAFFT, MUSCLE, Clustal Omega, T-Coffee, and ProbCons) on the BAliBASE 3.0 benchmark dataset to see how they stack up in terms of accuracy and performance.

## What This Does

This project runs several MSA tools on the same benchmark sequences and measures:
- Alignment accuracy (SP and TC scores)
- Runtime and memory usage
- Performance across different types of sequence families

Everything runs in Docker so you don't have to manually install the MSA tools.

## Quick Start

```bash
# Clone the repo
git clone https://github.com/ibrqures-uf/compass.git
cd compass

# Build and run
docker build -t msa-benchmark .
docker run -v "${PWD}:/app" -w /app msa-benchmark python3 main.py
```

The full benchmark takes a while to run (386 sequence families × 5 tools = 1,930 alignments). To test it out first:

```bash
docker run -e BENCH_LIMIT=5 -v "${PWD}:/app" -w /app msa-benchmark python3 main.py
```

## What You Need

- Docker installed on your machine
- At least 4GB RAM available
- About 2GB disk space for Docker image and results

## How It Works

The pipeline has four stages:

1. **Data Setup** - Downloads BAliBASE benchmark sequences (or uses existing data)
2. **Run Alignments** - Runs all five MSA tools on each test case
3. **Score Results** - Compares alignments to BAliBASE reference using SP/TC metrics
4. **Visualize** - Generates comparison plots and summary statistics

## Output

Results go in the `results/` directory:
- `benchmark_results.csv` - All the raw data
- `figures/` - Comparison plots showing which tools perform best
- `alignments/` - The actual alignment files generated

## The MSA Tools We Test

- **MAFFT** - Fast progressive alignment with iterative refinement
- **MUSCLE** - Multiple sequence comparison by log-expectation
- **Clustal Omega** - Uses HMM profile-profile techniques
- **T-Coffee** - Consistency-based progressive alignment
- **ProbCons** - Probabilistic consistency-based alignment

All tools are installed automatically via Ubuntu packages in the Docker container.

## What We Found

Check out our report for the full analysis, but the short version:
- MUSCLE is fastest and uses least memory
- Clustal Omega has highest mean SP and TC scores
- T-Coffee is slowest but produces decent alignments
- Different tools excel at different types of sequences

## Project Structure

```
src/
├── main.py              # Main pipeline controller
├── api/
│   └── data_fetcher.py  # Downloads and manages BAliBASE data
├── tools/
│   └── msa_runner.py    # Runs MSA tools and monitors resources
├── scoring/
│   └── baliscore.py     # Calculates SP/TC scores
└── analysis/
    ├── aggregator.py    # Combines results
    └── visualizer.py    # Makes plots
```

## Notes

- This is for benchmarking existing tools, not developing new algorithms
- We use standard BAliBASE scoring (SP and TC metrics)
- The Docker setup ensures consistent tool versions across different machines
- You'll need to download BAliBASE 3.0 separately if you want the actual data (we create placeholder structure otherwise)

## CAP5510 Final Project

This was our final project for CAP5510 Bioinformatics. The goal was to create a reproducible benchmarking pipeline and compare MSA tool performance.

**Team**: Aelly Alwardi, Daen Khan Patan, Ibraheem Qureshi
