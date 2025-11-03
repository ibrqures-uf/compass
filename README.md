# MSA Benchmarking Suite

A containerized benchmarking system for Multiple Sequence Alignment (MSA) tools using the BAliBASE benchmark dataset.

## 🚀 Quick Start

```bash
# Build the Docker image
docker build -t msa-benchmark .

# Run the benchmark
docker run -v "${PWD}:/app" -w /app msa-benchmark python3 main.py
```

## 📋 Features

- **Automated Tool Installation**: MAFFT, MUSCLE, and Clustal Omega are automatically installed
- **BAliBASE Integration**: Automatic download and processing of BAliBASE benchmark datasets
- **Multiple Format Support**: Handles FASTA, MSF, and RSF alignment formats
- **Comprehensive Scoring**: Calculates both SP (Sum-of-Pairs) and TC (Total Column) scores
- **Result Visualization**: Generates performance comparison plots and summary statistics
- **Docker-Based**: Runs entirely in container with no host dependencies

## 🛠️ Prerequisites

- Docker Desktop (Windows/macOS) or Docker Engine (Linux)
- At least 4GB of available RAM
- 2GB of free disk space

## 💻 Installation

1. Clone this repository:
   ```bash
   git clone <repository-url>
   cd bio_informatics
   ```

2. Build the Docker image:
   ```bash
   docker build -t msa-benchmark .
   ```

## 🚀 Usage

### Basic Run
```bash
docker run -v "${PWD}:/app" -w /app msa-benchmark python3 main.py
```

### With Limited Dataset (Testing)
```bash
docker run -e BENCH_LIMIT=5 -v "${PWD}:/app" -w /app msa-benchmark python3 main.py
```

### With Resource Limits
```bash
docker run --memory=4g --cpus=2 -v "${PWD}:/app" -w /app msa-benchmark python3 main.py
```

## 📊 Output

The benchmark generates several outputs in the `results/` directory:

- `results/benchmark_results.csv`: Raw benchmark data
- `results/alignments/`: Generated MSA files
- `results/figures/`:
  - `accuracy_comparison.png`: SP/TC score comparison
  - `efficiency_comparison.png`: Runtime/memory usage
  - `performance_by_refset.png`: Performance across reference sets

## 📈 Scoring Metrics

- **SP Score (Sum-of-Pairs)**: Measures alignment accuracy by comparing aligned residue pairs
- **TC Score (Total Column)**: Measures the fraction of correctly aligned columns
- **Runtime**: Execution time in seconds
- **Memory Usage**: Peak memory usage in MB

## 🔧 Configuration

### Environment Variables

- `BENCH_LIMIT`: Limit the number of sequences to process (e.g., `5` for testing)
- `PYTHONPATH`: Automatically set by Docker to `/app`

### Resource Recommendations

- **Minimal**: 2GB RAM, 1 CPU
- **Recommended**: 4GB RAM, 2 CPUs
- **Full Dataset**: 8GB RAM, 4 CPUs

## 📝 Supported MSA Tools

| Tool | Version | Status |
|------|---------|---------|
| MAFFT | Latest | ✅ Included |
| MUSCLE | Latest | ✅ Included |
| Clustal Omega | Latest | ✅ Included |
| T-Coffee | - | ⚠️ Optional |
| ProbCons | - | ⚠️ Optional |

## 🤝 Contributing

Contributions are welcome! Please feel free to submit pull requests.

## 📜 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 🙏 Acknowledgments

- BAliBASE dataset providers
- Developers of MAFFT, MUSCLE, and Clustal Omega
- Python Bio community

## 📞 Support

For issues and questions:
1. Create an issue in the repository
2. Include detailed reproduction steps
3. Attach relevant error messages and logs
