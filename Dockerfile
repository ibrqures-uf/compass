# Use Ubuntu as base image for better tool compatibility
FROM ubuntu:22.04

# Avoid prompts from apt
ENV DEBIAN_FRONTEND=noninteractive

# Set working directory
WORKDIR /app


# Install system dependencies (including 'time' for memory measurement)
RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip \
    wget \
    curl \
    tar \
    unzip \
    build-essential \
    zlib1g-dev \
    time \
    && rm -rf /var/lib/apt/lists/*

# Install MSA tools
RUN apt-get update && apt-get install -y \
    mafft \
    muscle \
    clustalo \
    t-coffee \
    probcons \
    && rm -rf /var/lib/apt/lists/*

# Verify tool installations
RUN which mafft && \
    which muscle && \
    which clustalo && \
    which t_coffee && \
    which probcons

# Copy project files
COPY . /app/

# Install Python dependencies
RUN pip3 install --no-cache-dir -r requirements.txt

# Create necessary directories
RUN mkdir -p data/balibase results/alignments results/figures

# Set environment variables
ENV PYTHONPATH=/app
ENV PATH="/app/tools:${PATH}"

# Default command
CMD ["python3", "main.py"]