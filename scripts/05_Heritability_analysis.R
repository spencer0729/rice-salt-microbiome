# ============================================================
# Script: 05_Heritability_analysis.R
# Purpose: Estimate broad-sense heritability of ASVs using linear mixed models.
# ============================================================

library(lme4)
library(dplyr)

# =========================
# Input data
# =========================
# The input table should contain:
# line: rice genotype
# loc : field/site
# rep : biological replicate within site
# ASV columns from the 4th column onward

phe <- read.table(
  "rhizos_individual.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

phe$line <- as.factor(phe$line)
phe$loc  <- as.factor(phe$loc)
phe$rep  <- as.factor(phe$rep)

asv_cols <- colnames(phe)[4:ncol(phe)]

# =========================
# Function for H2 calculation
# =========================

calc_H2 <- function(trait_name, data) {
  
  model <- lmer(
    as.formula(
      paste0(
        "`", trait_name, "` ~ ",
        "(1|line) + (1|loc) + (1|rep%in%loc) + (1|line:loc)"
      )
    ),
    data = data,
    control = lmerControl(
      check.nobs.vs.rankZ = "warning",
      check.nobs.vs.nlev  = "warning",
      check.nobs.vs.nRE   = "warning",
      check.nlev.gtreq.5  = "warning",
      check.nlev.gtr.1    = "warning"
    )
  )
  
  vc <- as.data.frame(VarCorr(model))
  
  var_line <- vc$vcov[vc$grp == "line"]
  var_line_loc <- vc$vcov[vc$grp == "line:loc"]
  var_residual <- vc$vcov[vc$grp == "Residual"]
  
  H2 <- var_line / (
    var_line +
      var_line_loc / 2 +
      var_residual / 4
  )
  
  data.frame(
    trait = trait_name,
    var_line = var_line,
    var_line_loc = var_line_loc,
    var_residual = var_residual,
    H2 = H2
  )
}

# =========================
# Batch calculation
# =========================

H2_result <- lapply(
  asv_cols,
  calc_H2,
  data = phe
) %>%
  bind_rows()

# =========================
# Output
# =========================

write.table(
  H2_result,
  "ASV_heritability_H2_result.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Heritable ASVs
heritable_ASVs <- H2_result %>%
  filter(H2 > 0.5)

write.table(
  heritable_ASVs,
  "heritable_ASVs_H2_gt_0.5.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)