#!/usr/bin/env Rscript

# ============================================================
# Script: 06_Env_STI_microbiome_analysis.R
# Purpose: Perform Mantel tests, Spearman correlation and variation partitioning for environment, STI and microbiome data.
# ============================================================

library(vegan)
library(psych)
library(reshape2)
library(ggplot2)

# =========================
# 1. Read data
# =========================

env <- read.delim("env_genji.txt", sep = "\t", row.names = 1)
STI <- read.delim("STI_genji.txt", sep = "\t", row.names = 1)
spe <- read.delim("asvtab_chouping_genji.txt", sep = "\t", row.names = 1)

spe <- t(spe)
spe <- spe[rownames(env), ]

# =========================
# 2. Variation partitioning
# =========================

spe_hel <- decostand(spe, method = "hellinger")

rda_vp <- varpart(spe_hel, env, STI)
pdf("variation_partitioning_env_STI.pdf", width = 4, height = 4)
plot(rda_vp, digits = 2, Xnames = c("Environment", "STI"),
     bg = c("#4DBBD5", "#E64B35"))
dev.off()

env_test <- anova.cca(rda(spe_hel, env), permutations = 999)
sti_test <- anova.cca(rda(spe_hel, STI), permutations = 999)

write.table(env_test, "RDA_env_permutation_test.txt", sep = "\t", quote = FALSE)
write.table(sti_test, "RDA_STI_permutation_test.txt", sep = "\t", quote = FALSE)

# =========================
# 3. Mantel test
# =========================

dist_spe <- vegdist(spe_hel, method = "bray")
dist_STI <- dist(STI$STI, method = "euclidean")

env_scaled <- scale(env, center = TRUE, scale = TRUE)
dist_env <- dist(env_scaled, method = "euclidean")

mantel_STI <- mantel(dist_spe, dist_STI, method = "spearman",
                     permutations = 9999, na.rm = TRUE)

mantel_env <- mantel(dist_spe, dist_env, method = "spearman",
                     permutations = 9999, na.rm = TRUE)

mantel_result <- data.frame(
  variable = c("STI", "Environment"),
  mantel_r = c(mantel_STI$statistic, mantel_env$statistic),
  p_value = c(mantel_STI$signif, mantel_env$signif)
)

write.table(mantel_result, "mantel_result.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# =========================
# 4. Spearman correlation
# =========================

dat <- cbind(STI, spe)
dat <- decostand(dat, method = "normalize")

STI_norm <- dat[, "STI", drop = FALSE]
spe_norm <- dat[, colnames(dat) != "STI"]

spearman <- corr.test(STI_norm, spe_norm,
                      method = "spearman",
                      adjust = "fdr")

r <- data.frame(spearman$r)
p <- data.frame(spearman$p)

r$env <- rownames(r)
p$env <- rownames(p)

r <- melt(r, id = "env")
p <- melt(p, id = "env")

spearman_result <- cbind(r, p$value)
colnames(spearman_result) <- c("env", "ASV",
                               "spearman_correlation", "p_value")

spearman_result$sig <- ""
spearman_result$sig[spearman_result$p_value < 0.05] <- "*"
spearman_result$sig[spearman_result$p_value < 0.01] <- "**"
spearman_result$sig[spearman_result$p_value < 0.001] <- "***"

write.table(spearman_result, "spearman_STI_ASV_result.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

