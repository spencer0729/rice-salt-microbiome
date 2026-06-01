#!/usr/bin/env bash
# Configuration file for 01_QIIME2_DADA2_pipeline.sh

RAW_DIR="data/raw_fastq"
RESULTS_DIR="results/qiime2"
METADATA="metadata.tsv"
CLASSIFIER="silva-138-99-nb-classifier.qza"

THREADS=0
TRIM_LEFT_F=0
TRIM_LEFT_R=0
TRUNC_LEN_F=250
TRUNC_LEN_R=250

EXCLUDE_TAXA="mitochondria,chloroplast,eukaryota"

mkdir -p "${RESULTS_DIR}"
