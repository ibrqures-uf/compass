FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

WORKDIR /app

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

RUN apt-get update && apt-get install -y \
    mafft \
    muscle \
    clustalo \
    t-coffee \
    probcons \
    && rm -rf /var/lib/apt/lists/*

RUN which mafft && \
    which muscle && \
    which clustalo && \
    which t_coffee && \
    which probcons

COPY . /app/

RUN pip3 install --no-cache-dir -r requirements.txt

RUN mkdir -p data/balibase results/alignments results/figures

ENV PYTHONPATH=/app
ENV PATH="/app/tools:${PATH}"

CMD ["python3", "main.py"]
