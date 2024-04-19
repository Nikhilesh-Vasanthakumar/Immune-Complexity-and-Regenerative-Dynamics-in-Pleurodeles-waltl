# Requirements for Conda Environment

This document describes the requirements for setting up a Conda environment with the necessary packages for our bioinformatics pipeline.

## Environment Setup

To create a new Conda environment and install the required packages, follow these steps:

### 1. Create a new Conda environment

conda create -n bioinfo_env python=3.9

### 2. Activate the environment

conda activate bioinfo_env

### 3. Install packages

conda install -c bioconda scanpy=1.9.6 anndata=0.9.2 scikit-learn=1.3.1 numpy=1.23.5
conda install -c conda-forge harmony-pytorch=0.0.9
pip install samap

### Define NCBI BLAST version.
ncbi_blast_version='2.9.0'

### Download NCBI BLAST tarball.
wget "ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${ncbi_blast_version}/ncbi-blast-${ncbi_blast_version}+-x64-linux.tar.gz"

### Extract NCBI BLAST binaries in current conda environment bin directory.
tar -xzvf "ncbi-blast-${ncbi_blast_version}+-x64-linux.tar.gz" \
    -C "${CONDA_PREFIX}/bin/" \
    --strip-components=2 \
    "ncbi-blast-${ncbi_blast_version}+/bin/"
