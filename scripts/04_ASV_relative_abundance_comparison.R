# ============================================================
# Script: 04_ASV_relative_abundance_comparison.R
# Purpose: Compare relative abundance of selected ASVs among groups using ANOVA and LSD tests.
# ============================================================

library(ggplot2)
library(dplyr)
library(agricolae)
library(patchwork)

# Input files
asv_files <- c(
  "asv_22.csv",
  "asv_54.csv",
  "asv_60.csv",
  "asv_96.csv",
  "asv_107.csv",
  "asv_114.csv",
  "asv_155.csv"
)

asv_names <- c(
  "ASV22",
  "ASV54",
  "ASV60",
  "ASV96",
  "ASV107",
  "ASV114",
  "ASV155"
)

group_order <- c("NS_R", "SS_R", "NS_E", "SS_E")

group_cols <- c(
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
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 7)
  )

plot_asv <- function(file, title) {
  
  df <- read.csv(file, stringsAsFactors = FALSE, row.names = 1)
  df$group <- factor(df$group, levels = group_order)
  
  model <- aov(relative.abundance ~ group, data = df)
  lsd <- LSD.test(model, "group", p.adj = "none")
  
  stat <- lsd$groups
  df$stat <- stat[as.character(df$group), "groups"]
  
  y_pos <- df %>%
    group_by(group) %>%
    summarise(y = max(relative.abundance, na.rm = TRUE), .groups = "drop")
  
  y_range <- max(df$relative.abundance, na.rm = TRUE) -
    min(df$relative.abundance, na.rm = TRUE)
  
  y_pos$y <- y_pos$y + y_range * 0.08
  df <- left_join(df, y_pos, by = "group")
  
  p <- ggplot(df, aes(group, relative.abundance, color = group)) +
    geom_boxplot(width = 0.5, fill = "transparent", outlier.shape = NA, linewidth = 0.4) +
    geom_jitter(width = 0.17, size = 0.7, alpha = 0.7) +
    geom_text(aes(y = y, label = stat), size = 7 / .pt) +
    scale_x_discrete(limits = group_order) +
    scale_color_manual(values = group_cols) +
    labs(x = NULL, y = "Relative Abundance(%)", title = title) +
    main_theme
  
  return(p)
}

plot_list <- mapply(
  plot_asv,
  asv_files,
  asv_names,
  SIMPLIFY = FALSE
)

p_all <- wrap_plots(plot_list, ncol = 4)

ggsave(
  "ASV_relative_abundance_comparison.pdf",
  p_all,
  width = 8,
  height = 4,
  device = cairo_pdf
)