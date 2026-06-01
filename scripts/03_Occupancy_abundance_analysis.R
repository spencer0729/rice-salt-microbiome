# ============================================================
# Script: 03_Occupancy_abundance_analysis.R
# Purpose: Calculate occupancy and mean relative abundance of ASVs and generate occupancy-abundance plot.
# ============================================================

library(tidyverse)
library(vegan)

# Read data
otu <- read.table("otu_VE.txt", header = TRUE, row.names = 1, sep = "\t")

# Calculate occupancy and mean relative abundance
otu_PA <- ifelse(otu > 0, 1, 0)
otu_occ <- rowSums(otu_PA) / ncol(otu)
otu_rel <- apply(decostand(otu, method = "total", MARGIN = 2), 1, mean)

occ_abun <- data.frame(
  ASV = rownames(otu),
  occupancy = otu_occ,
  mean_relative_abundance = otu_rel
)

write.table(
  occ_abun,
  "occupancy_abundance_result.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Read annotated table
# The table should contain: ASV, occupancy, mean_relative_abundance, fill
occ_abun <- read.table(
  "occupancy_abundance_annotated.txt",
  header = TRUE,
  sep = "\t"
)

occ_abun$fill <- factor(
  occ_abun$fill,
  levels = c("Core ASVs", "Non-core ASVs")
)

# Plot
p <- ggplot(
  occ_abun,
  aes(
    x = log10(mean_relative_abundance),
    y = occupancy,
    color = fill,
    shape = fill
  )
) +
  geom_point(size = 1, alpha = 0.8) +
  scale_color_manual(
    values = c(
      "Core ASVs" = "#F98B60",
      "Non-core ASVs" = "#CDCDC1"
    ),
    name = NULL
  ) +
  scale_shape_manual(
    values = c(
      "Core ASVs" = 19,
      "Non-core ASVs" = 1
    ),
    name = NULL
  ) +
  labs(
    x = paste0(
      "log10(mean relative abundance per ASV)\n(n = ",
      nrow(occ_abun),
      " ASVs)"
    ),
    y = paste0(
      "Occupancy (n = ",
      ncol(otu),
      ")"
    )
  ) +
  theme_classic(base_size = 7) +
  theme(
    text = element_text(family = "Arial", size = 7),
    axis.text = element_text(size = 7, colour = "black"),
    axis.title = element_text(size = 7, colour = "black"),
    legend.text = element_text(size = 7),
    legend.position = "right"
  )

ggsave(
  "Fig_occupancy_abundance.pdf",
  p,
  width = 4,
  height = 2.5,
  device = cairo_pdf
)