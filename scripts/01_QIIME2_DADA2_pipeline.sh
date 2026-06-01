#!/usr/bin/env bash
# ============================================================
# Script: 01_QIIME2_DADA2_pipeline.sh
# Purpose: QIIME2 processing pipeline for paired-end 16S rRNA sequencing reads.
# ============================================================

set -euo pipefail

source "$(dirname "$0")/../config/config.sh"

MANIFEST="${RESULTS_DIR}/manifest.tsv"

if [[ ! -f "${MANIFEST}" ]]; then
  echo "Manifest not found: ${MANIFEST}"
  echo "Run: Rscript scripts/00_create_manifest.R data/raw_fastq results/manifest.tsv"
  exit 1
fi

qiime tools import   --type 'SampleData[PairedEndSequencesWithQuality]'   --input-path "${MANIFEST}"   --output-path "${RESULTS_DIR}/demux.qza"   --input-format PairedEndFastqManifestPhred33V2

qiime demux summarize   --i-data "${RESULTS_DIR}/demux.qza"   --o-visualization "${RESULTS_DIR}/demux.qzv"

qiime dada2 denoise-paired   --i-demultiplexed-seqs "${RESULTS_DIR}/demux.qza"   --p-n-threads "${THREADS}"   --p-trim-left-f "${TRIM_LEFT_F}"   --p-trim-left-r "${TRIM_LEFT_R}"   --p-trunc-len-f "${TRUNC_LEN_F}"   --p-trunc-len-r "${TRUNC_LEN_R}"   --o-table "${RESULTS_DIR}/table.qza"   --o-representative-sequences "${RESULTS_DIR}/rep-seqs.qza"   --o-denoising-stats "${RESULTS_DIR}/denoising-stats.qza"

qiime metadata tabulate   --m-input-file "${RESULTS_DIR}/denoising-stats.qza"   --o-visualization "${RESULTS_DIR}/denoising-stats.qzv"

qiime feature-table tabulate-seqs   --i-data "${RESULTS_DIR}/rep-seqs.qza"   --o-visualization "${RESULTS_DIR}/rep-seqs.qzv"

qiime feature-table summarize   --i-table "${RESULTS_DIR}/table.qza"   --o-visualization "${RESULTS_DIR}/table.qzv"   --m-sample-metadata-file "${METADATA}"

qiime feature-classifier classify-sklearn   --i-classifier "${CLASSIFIER}"   --i-reads "${RESULTS_DIR}/rep-seqs.qza"   --o-classification "${RESULTS_DIR}/taxonomy.qza"

qiime metadata tabulate   --m-input-file "${RESULTS_DIR}/taxonomy.qza"   --o-visualization "${RESULTS_DIR}/taxonomy.qzv"

qiime taxa filter-table   --i-table "${RESULTS_DIR}/table.qza"   --i-taxonomy "${RESULTS_DIR}/taxonomy.qza"   --p-exclude "${EXCLUDE_TAXA}"   --o-filtered-table "${RESULTS_DIR}/table-filtered.qza"

qiime feature-table summarize   --i-table "${RESULTS_DIR}/table-filtered.qza"   --o-visualization "${RESULTS_DIR}/table-filtered.qzv"   --m-sample-metadata-file "${METADATA}"

qiime tools export   --input-path "${RESULTS_DIR}/table-filtered.qza"   --output-path "${RESULTS_DIR}/export"

biom convert   -i "${RESULTS_DIR}/export/feature-table.biom"   -o "${RESULTS_DIR}/export/otu_table.tsv"   --to-tsv

sed -i '1d' "${RESULTS_DIR}/export/otu_table.tsv"
sed -i 's/#OTU ID/ASV/' "${RESULTS_DIR}/export/otu_table.tsv"

qiime tools export   --input-path "${RESULTS_DIR}/taxonomy.qza"   --output-path "${RESULTS_DIR}/export"

qiime tools export   --input-path "${RESULTS_DIR}/rep-seqs.qza"   --output-path "${RESULTS_DIR}/export"

qiime phylogeny align-to-tree-mafft-fasttree   --i-sequences "${RESULTS_DIR}/rep-seqs.qza"   --o-alignment "${RESULTS_DIR}/aligned-rep-seqs.qza"   --o-masked-alignment "${RESULTS_DIR}/masked-aligned-rep-seqs.qza"   --o-tree "${RESULTS_DIR}/unrooted-tree.qza"   --o-rooted-tree "${RESULTS_DIR}/rooted-tree.qza"