#!/usr/bin/env Rscript

# ============================================================
# Script: 02_Microbiome_diversity.R
# Purpose: Calculate alpha diversity, Bray-Curtis distance, PCoA coordinates and PERMANOVA results.
# ============================================================

# 01_calculate_alpha_beta_diversity.R
# Purpose: Calculate alpha diversity indices, Bray-Curtis distance, PCoA coordinates,
#          and PERMANOVA results from an ASV table.
# Input files:
#   - asvtab.csv: ASV table, rows = ASVs/features, columns = samples
#   - group.csv: sample metadata. Required columns: SampleID, group
# Output folders:
#   - data_result/
#   - alpha_figs/
#   - beta_figs/

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(ggplot2)
})

set.seed(123)

# -------------------------
# 1. Create output folders
# -------------------------
dir.create("data_result", showWarnings = FALSE, recursive = TRUE)
dir.create("alpha_figs", showWarnings = FALSE, recursive = TRUE)
dir.create("beta_figs", showWarnings = FALSE, recursive = TRUE)

# -------------------------
# 2. Read input data
# -------------------------
metadata <- read.csv("group.csv", header = TRUE, stringsAsFactors = FALSE)
asv_table <- read.csv("asvtab.csv", header = TRUE, row.names = 1, check.names = FALSE) %>%
  t() %>%
  as.data.frame()

stopifnot("SampleID" %in% colnames(metadata))
stopifnot("group" %in% colnames(metadata))

# Keep only samples shared by ASV table and metadata
shared_samples <- intersect(rownames(asv_table), metadata$SampleID)
asv_table <- asv_table[shared_samples, , drop = FALSE]
metadata <- metadata %>%
  filter(SampleID %in% shared_samples) %>%
  arrange(match(SampleID, rownames(asv_table)))
rownames(metadata) <- metadata$SampleID

# -------------------------
# 3. Rarefy ASV table
# -------------------------
min_depth <- min(rowSums(asv_table))
asv_rare <- rrarefy(asv_table, sample = min_depth)
write.csv(asv_rare, "data_result/asv_table_rarefied.csv")

# Optional: rarefaction curve
pdf("alpha_figs/rarefaction_curve.pdf", width = 7, height = 5)
rarecurve(asv_rare, step = 2000, col = rainbow(nrow(asv_rare)), label = FALSE)
dev.off()

# -------------------------
# 4. Alpha diversity
# -------------------------
alpha_index <- data.frame(
  SampleID = rownames(asv_rare),
  observed = specnumber(asv_rare),
  chao1 = estimateR(t(asv_rare))["S.chao1", ],
  ace = estimateR(t(asv_rare))["S.ACE", ],
  shannon = diversity(asv_rare, index = "shannon"),
  simpson = diversity(asv_rare, index = "simpson"),
  group = metadata[rownames(asv_rare), "group"],
  row.names = rownames(asv_rare)
)
write.csv(alpha_index, "data_result/alpha_index.csv")

# -------------------------
# 5. Bray-Curtis distance and PCoA
# -------------------------
bray_dist <- vegdist(asv_rare, method = "bray")
write.table(as.matrix(bray_dist), "data_result/bray_curtis_distance.tsv",
            sep = "\t", quote = FALSE, col.names = NA)

pcoa <- cmdscale(bray_dist, k = min(nrow(asv_rare) - 1, 10), eig = TRUE)
var_explained <- 100 * pcoa$eig / sum(pcoa$eig[pcoa$eig > 0])

pcoa_points <- as.data.frame(pcoa$points)
colnames(pcoa_points) <- paste0("PCoA", seq_len(ncol(pcoa_points)))
pcoa_points <- cbind(SampleID = rownames(pcoa_points), pcoa_points, metadata[rownames(pcoa_points), , drop = FALSE])

write.csv(pcoa_points, "data_result/PCoA_points.csv", row.names = FALSE)
write.csv(data.frame(axis = paste0("PCoA", seq_along(var_explained)),
                     variance_explained = var_explained),
          "data_result/PCoA_variance_explained.csv", row.names = FALSE)

# Scree plot
scree_df <- data.frame(axis = seq_along(var_explained), variance_explained = var_explained)
p_scree <- ggplot(scree_df, aes(x = axis, y = variance_explained)) +
  geom_point(size = 1.8) +
  geom_line(linewidth = 0.4) +
  theme_classic(base_size = 7) +
  labs(x = "PCoA axis", y = "Variance explained (%)")
ggsave("beta_figs/PCoA_scree_plot.pdf", p_scree, width = 3.5, height = 2.5)

# -------------------------
# 6. PERMANOVA
# -------------------------
permanova_all <- adonis2(bray_dist ~ group, data = metadata, permutations = 999)
write.table(as.data.frame(permanova_all), "data_result/PERMANOVA_all_groups.tsv",
            sep = "\t", quote = FALSE, col.names = NA)

# Pairwise PERMANOVA
pairwise_permanova <- function(asv_mat, meta, group_col = "group") {
  groups <- unique(meta[[group_col]])
  out <- list()
  k <- 1
  for (i in seq_len(length(groups) - 1)) {
    for (j in (i + 1):length(groups)) {
      selected_groups <- c(groups[i], groups[j])
      sub_meta <- meta[meta[[group_col]] %in% selected_groups, , drop = FALSE]
      sub_asv <- asv_mat[sub_meta$SampleID, , drop = FALSE]
      sub_dist <- vegdist(sub_asv, method = "bray")
      fit <- adonis2(sub_dist ~ group, data = sub_meta, permutations = 999)
      out[[k]] <- data.frame(
        comparison = paste(selected_groups, collapse = " vs "),
        R2 = fit$R2[1],
        F = fit$F[1],
        P = fit$`Pr(>F)`[1]
      )
      k <- k + 1
    }
  }
  bind_rows(out) %>% mutate(P_adj_BH = p.adjust(P, method = "BH"))
}

pairwise_result <- pairwise_permanova(asv_rare, metadata)
write.table(pairwise_result, "data_result/PERMANOVA_pairwise.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# 02_plot_shannon_violin.R
# Purpose: Draw Shannon diversity violin plot.
# Input file:
#   - alpha_index.csv or data_result/alpha_index.csv
# Required columns:
#   - group
#   - shannon

suppressPackageStartupMessages({
  library(tidyverse)
  library(rstatix)
  library(ggpubr)
  library(ggplot2)
})

# -------------------------
# 1. Parameters
# -------------------------
input_file <- if (file.exists("data_result/alpha_index.csv")) {
  "data_result/alpha_index.csv"
} else {
  "alpha_index.csv"
}

output_file <- "Shannon_violin.pdf"
group_levels <- c("NS_R", "SS_R", "NS_E", "SS_E")
comparisons <- list(c("NS_R", "SS_R"), c("NS_E", "SS_E"))

cols_group <- c(
  "NS_R" = "#66a61e",
  "SS_R" = "#8e7cc3",
  "NS_E" = "#e67e73",
  "SS_E" = "#1abc9c"
)

main_theme <- theme_classic(base_size = 7) +
  theme(
    text = element_text(family = "Arial", size = 7),
    axis.text = element_text(size = 7, colour = "black"),
    axis.title = element_text(size = 7, colour = "black"),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    legend.position = "none"
  )

# -------------------------
# 2. Read data
# -------------------------
df <- read.csv(input_file, stringsAsFactors = FALSE, row.names = 1)
stopifnot(all(c("group", "shannon") %in% colnames(df)))

df <- df %>%
  filter(group %in% group_levels) %>%
  mutate(group = factor(group, levels = group_levels))

# -------------------------
# 3. Welch t-test
# -------------------------
p_dat <- df %>%
  t_test(shannon ~ group, comparisons = comparisons, var.equal = FALSE) %>%
  add_significance("p")

y_top <- max(df$shannon, na.rm = TRUE)
p_dat$y.position <- c(y_top + 0.12, y_top + 0.30)

# -------------------------
# 4. Plot
# -------------------------
p <- ggplot(df, aes(x = group, y = shannon, fill = group)) +
  geom_violin(width = 0.9, trim = FALSE, color = NA, alpha = 0.35) +
  geom_jitter(aes(color = group), width = 0.12, alpha = 0.55, size = 1, show.legend = FALSE) +
  stat_summary(fun = mean, geom = "point", size = 2, color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15, linewidth = 0.5, color = "black") +
  stat_pvalue_manual(
    p_dat,
    label = "p.signif",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    size = 7 / .pt,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = cols_group) +
  scale_color_manual(values = cols_group) +
  scale_y_continuous(breaks = seq(1, 6, 1)) +
  coord_cartesian(ylim = c(min(df$shannon, na.rm = TRUE) - 0.2, y_top + 0.45)) +
  labs(x = NULL, y = "Shannon") +
  main_theme

ggsave(output_file, p, width = 3.5, height = 2.5, device = cairo_pdf)

# 03_plot_pcoa_sti.R
# Purpose: Draw Bray-Curtis PCoA plot colored by STI.
# Input options:
#   Option A: asvtab.csv + design.txt
#   Option B: data_result/bray_curtis_distance.tsv + design.txt
# Required metadata columns in design.txt:
#   - group: NS_R, SS_R, NS_E, SS_E
#   - STI: salt-tolerance index
# Optional metadata columns:
#   - compartment: R/E or rhizosphere/endosphere
#   - site: NS/SS or normal/saline field

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(ggplot2)
})

set.seed(123)

# -------------------------
# 1. Parameters
# -------------------------
metadata_file <- "design.txt"
asv_file <- "asvtab.csv"
distance_file <- "data_result/bray_curtis_distance.tsv"
output_file <- "PCoA_STI.pdf"
group_levels <- c("NS_R", "SS_R", "NS_E", "SS_E")

main_theme <- theme(
  panel.background = element_blank(),
  panel.grid = element_blank(),
  axis.line.x = element_line(linewidth = 0.5, colour = "black"),
  axis.line.y = element_line(linewidth = 0.5, colour = "black"),
  axis.ticks = element_line(color = "black", linewidth = 0.5),
  axis.text = element_text(color = "black", size = 7),
  axis.title = element_text(color = "black", size = 7),
  legend.position = "right",
  legend.background = element_blank(),
  legend.key = element_blank(),
  legend.text = element_text(size = 7),
  legend.title = element_text(size = 7),
  text = element_text(family = "Arial", size = 7)
)

# -------------------------
# 2. Read metadata
# -------------------------
design <- read.table(metadata_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
stopifnot("group" %in% colnames(design))
stopifnot("STI" %in% colnames(design))

design <- design %>%
  rownames_to_column("SampleID") %>%
  filter(group %in% group_levels) %>%
  mutate(group = factor(group, levels = group_levels)) %>%
  column_to_rownames("SampleID")

# -------------------------
# 3. Calculate or read Bray-Curtis distance
# -------------------------
if (file.exists(distance_file)) {
  bray <- read.delim(distance_file, row.names = 1, check.names = FALSE)
} else {
  asv_table <- read.csv(asv_file, header = TRUE, row.names = 1, check.names = FALSE) %>%
    t() %>%
    as.data.frame()
  shared_samples <- intersect(rownames(asv_table), rownames(design))
  asv_table <- asv_table[shared_samples, , drop = FALSE]
  asv_rare <- rrarefy(asv_table, sample = min(rowSums(asv_table)))
  bray <- as.matrix(vegdist(asv_rare, method = "bray"))
}

# Match samples between distance matrix and metadata
shared_samples <- intersect(rownames(bray), rownames(design))
bray <- bray[shared_samples, shared_samples]
design <- design[shared_samples, , drop = FALSE]

# -------------------------
# 4. PCoA
# -------------------------
pcoa <- cmdscale(as.dist(bray), k = 4, eig = TRUE)
positive_eig <- pcoa$eig[pcoa$eig > 0]
var1 <- 100 * pcoa$eig[1] / sum(positive_eig)
var3 <- 100 * pcoa$eig[3] / sum(positive_eig)

points <- as.data.frame(pcoa$points)
colnames(points) <- c("PCoA1", "PCoA2", "PCoA3", "PCoA4")
points <- cbind(points, design[rownames(points), , drop = FALSE])

# Shape mapping: use group by default. This matches the published-style panel where
# rhizosphere/endosphere and field are encoded together as group.
points$group <- factor(points$group, levels = group_levels)

# -------------------------
# 5. Plot PCoA1 vs PCoA3 colored by STI
# -------------------------
p <- ggplot(points, aes(PCoA1, PCoA3, color = STI, shape = group)) +
  geom_point(alpha = 0.8, size = 1.6, stroke = 0.5) +
  scale_shape_manual(values = c(NS_R = 1, SS_R = 2, NS_E = 16, SS_E = 17), name = "Group") +
  scale_color_viridis_c(option = "D", name = "STI") +
  labs(x = paste0("PCoA 1 (", round(var1, 2), "%)"),
       y = paste0("PCoA 3 (", round(var3, 2), "%)")) +
  main_theme

ggsave("PCoA_STI.pdf", p, width = 4, height = 2.5, device = cairo_pdf)
