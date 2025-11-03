#!/bin/bash

# Update package list
apt-get update

# Install required tools
apt-get install -y \
    mafft \
    muscle \
    clustalo \
    t-coffee \
    probcons

# Verify installations
echo "Verifying installations..."
echo "MAFFT: $(which mafft)"
echo "MUSCLE: $(which muscle)"
echo "Clustal Omega: $(which clustalo)"
echo "T-Coffee: $(which t_coffee)"
echo "ProbCons: $(which probcons)"