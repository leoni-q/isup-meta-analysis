###############################################################################
# Evolution of ISUP Grade Group 1 Detection Across Prostate Cancer Screening
# Strategies: meta-analysis per protocol.md (section 8, Statistical Analysis)
#
# Models per protocol:
#   - Risk ratios: inverse-variance, DerSimonian-Laird tau, Hartung-Knapp CI,
#     0.5 continuity correction for zero cells, random effects.
#   - Proportions: study-strategy GLMMs (logit), fitted separately by strategy.
#   - Heterogeneity: Cochran's Q (p < 0.05 significant), plus tau^2 / I^2.
#   - Leave-one-study-out sensitivity analysis within every pooled group.
#   - Publication bias: NOT assessed (no analysis reaches k >= 10, the
#     protocol-specified threshold).
###############################################################################

# R syntax guide for readers coming from TS, Java, Python, or Rust:
#   <-              assignment (similar to = in most languages)
#   object$field    access a named field / data-frame column
#   object[[name]]  access a field whose name is stored in a variable
#   c(a, b)         create a vector
#   list(a = value) create a key-value list
#   data[rows, cols] select rows and columns (indices start at 1)
#   NA / is.na(x)   missing value / missing-value check
#   &, | / &&, ||   vectorized / single-condition logical operators
#
# This entrypoint runs from top to bottom. Each sourced module shares this
# script environment, so its behavior matches the previous single-file script.

ANALYSIS_DIR <- "/Users/lukaszkujawski/Code/ania/meta_analysis/analysis"
MODULE_DIR <- file.path(ANALYSIS_DIR, "R")

source(file.path(MODULE_DIR, "setup_and_data.R"), local = FALSE)
source(file.path(MODULE_DIR, "summary_helpers.R"), local = FALSE)
source(file.path(MODULE_DIR, "risk_ratio_helpers.R"), local = FALSE)
source(file.path(MODULE_DIR, "proportion_helpers.R"), local = FALSE)
source(file.path(MODULE_DIR, "part1_risk_ratios.R"), local = FALSE)
source(file.path(MODULE_DIR, "part2_proportions.R"), local = FALSE)
source(file.path(MODULE_DIR, "part3_sensitivity.R"), local = FALSE)
source(file.path(MODULE_DIR, "part4_finalize.R"), local = FALSE)
