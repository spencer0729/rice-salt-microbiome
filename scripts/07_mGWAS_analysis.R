# ============================================================
# Script: 07_mGWAS_analysis.R
# Purpose: Run/record GEMMA-based mGWAS workflow and generate Manhattan plots using CMplot.
# ============================================================

# ============================================================
# mGWAS_analysis.R
# Genome-wide association analysis and Manhattan plot
# ============================================================

library(CMplot)

# =========================
# 1. GWAS command lines
# =========================

# PCA analysis using PLINK:
# plink --bfile data --pca 20 --out data_pca

# Genomic relationship matrix using GEMMA:
# gemma -bfile data -gk 2 -p phenotype.txt -o phenotype_matrix

# Linear mixed model GWAS using GEMMA:
# gemma -bfile data \
#       -k output/phenotype_matrix.sXX.txt \
#       -lmm 1 \
#       -p phenotype.txt \
#       -c data_pca.eigenvec \
#       -o GWAS_result_phenotype

# Linkage disequilibrium analysis using LDBlockShow:
# LDBlockShow -InVCF data.vcf.gz \
#             -OutPut LD_region \
#             -Region chr:Start:end \
#             -SeleVar 2

# =========================
# 2. Read GWAS result
# =========================

gwas <- read.table(
  "GWAS_result.txt",
  header = TRUE,
  stringsAsFactors = FALSE
)

# The input file for CMplot should contain:
# SNP, Chromosome, Position, P.value

colnames(gwas)[1:4] <- c(
  "SNP",
  "Chromosome",
  "Position",
  "P.value"
)

# =========================
# 3. Significance threshold
# =========================

threshold <- 1.19e-08

sig_snp <- gwas$SNP[
  gwas$P.value < threshold
]

write.table(
  gwas[gwas$P.value < threshold, ],
  "significant_SNPs.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# =========================
# 4. Manhattan plot
# =========================

CMplot(
  gwas,
  plot.type = "m",
  LOG10 = TRUE,
  threshold = threshold,
  threshold.col = "blue",
  threshold.lty = 2,
  threshold.lwd = 1.5,
  amplify = TRUE,
  signal.cex = 1.2,
  signal.pch = 19,
  signal.col = "red",
  highlight = sig_snp,
  highlight.col = "red",
  highlight.cex = 1.3,
  highlight.pch = 19,
  chr.den.col = NULL,
  file = "pdf",
  dpi = 300,
  file.output = TRUE,
  verbose = TRUE,
  width = 12,
  height = 5
)

# =========================
# 5. Circular Manhattan plot
# =========================

CMplot(
  gwas,
  plot.type = "c",
  LOG10 = TRUE,
  outward = FALSE,
  threshold = threshold,
  threshold.col = "red",
  threshold.lty = 2,
  threshold.lwd = 1,
  chr.den.col = NULL,
  amplify = TRUE,
  signal.cex = 0.6,
  signal.pch = 19,
  cir.axis = TRUE,
  cir.axis.grid = TRUE,
  cir.chr = TRUE,
  file = "png",
  dpi = 600,
  file.output = TRUE,
  memo = "Circular_Manhattan"
)