# Rice microbiome and mGWAS analysis

This repository contains custom scripts used for microbiome diversity analysis,
ASV occupancy-abundance analysis, ASV heritability estimation, environmental association analysis,
and microbiome genome-wide association studies (mGWAS).

## Directory structure

```text
scripts/
├── 01_QIIME2_DADA2_pipeline.sh
├── 02_Microbiome_diversity.R
├── 03_Occupancy_abundance_analysis.R
├── 04_ASV_relative_abundance_comparison.R
├── 05_Heritability_analysis.R
├── 06_Env_STI_microbiome_analysis.R
└── 07_mGWAS_analysis.R

config/
└── config.sh
```

## Scripts

### 01_QIIME2_DADA2_pipeline.sh
QIIME2 workflow for paired-end 16S rRNA sequencing data, including importing reads, DADA2 denoising,
taxonomy assignment, filtering of unwanted taxa, feature-table export and phylogenetic tree construction.

### 02_Microbiome_diversity.R
Calculates alpha diversity indices, Bray-Curtis distance, PCoA coordinates, and PERMANOVA results.

### 03_Occupancy_abundance_analysis.R
Calculates ASV occupancy and mean relative abundance, and generates occupancy-abundance plots.

### 04_ASV_relative_abundance_comparison.R
Compares selected ASV relative abundance among NS_R, SS_R, NS_E and SS_E groups using ANOVA and LSD tests.

### 05_Heritability_analysis.R
Estimates broad-sense heritability (H²) of ASVs using linear mixed models.

### 06_Env_STI_microbiome_analysis.R
Performs variation partitioning, Mantel tests and Spearman correlation analyses among environmental variables,
salt-tolerance index (STI), and ASV profiles.

### 07_mGWAS_analysis.R
Provides the GEMMA/PLINK/LDBlockShow command-line workflow and uses CMplot for Manhattan and circular Manhattan plots.

## Notes

Input files are expected to be placed in the working directory or configured in `config/config.sh`.
Large sequencing files and intermediate QIIME2/GWAS output files are excluded from this repository.
