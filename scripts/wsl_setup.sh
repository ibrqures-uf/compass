#!/usr/bin/env bash
# WSL-side setup script: installs Miniconda, mamba, creates conda env with MSA tools,
# installs Python requirements and runs a small smoke-test benchmark.

set -euo pipefail

PROJECT_DIR="/mnt/c/Users/iqure/Desktop/bio_informatics"
echo "WSL setup: project dir is $PROJECT_DIR"

echo "Updating apt and installing prerequisites..."
sudo apt update
sudo apt install -y wget bzip2 ca-certificates curl git build-essential

MINICONDA_DIR="$HOME/miniconda"
if [ ! -d "$MINICONDA_DIR" ]; then
    echo "Installing Miniconda to $MINICONDA_DIR (non-interactive)..."
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh
    bash /tmp/miniconda.sh -b -p "$MINICONDA_DIR"
else
    echo "Miniconda already installed at $MINICONDA_DIR"
fi

export PATH="$MINICONDA_DIR/bin:$PATH"

echo "Initializing conda..."
eval "$(conda shell.bash hook)"

echo "Installing mamba into base environment..."
conda install -n base -c conda-forge mamba -y

echo "Creating conda environment 'msa_env' and installing MSA tools (bioconda)..."
# Use mamba for faster, more reliable installs
mamba create -n msa_env -c conda-forge -c bioconda mafft muscle clustalo t-coffee probcons -y

echo "Installing Python requirements into msa_env..."
conda run -n msa_env python -m pip install -r "$PROJECT_DIR/requirements.txt"

echo "Running a small smoke test (3 families) from the project..."
conda run -n msa_env python -m tools.benchmark_runner --data-dir "$PROJECT_DIR/data/balibase/bb3_release" --max-families 3 --out "$PROJECT_DIR/results/benchmark_results.csv"

echo "Smoke test finished. Results (if any) are at: $PROJECT_DIR/results/benchmark_results.csv"
echo "If you want to run the full benchmark:"
echo "  conda run -n msa_env python -m tools.benchmark_runner --data-dir $PROJECT_DIR/data/balibase/bb3_release --out $PROJECT_DIR/results/benchmark_results.csv"
