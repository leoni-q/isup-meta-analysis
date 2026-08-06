# Package setup, output paths, logging, and data preparation.

suppressPackageStartupMessages({
  library(meta)
  library(metafor)
})

settings.meta(digits = 3, digits.pval = 4)

OUTPUT_DIR <- file.path(ANALYSIS_DIR, "output")
RESULTS_LOG_PATH <- file.path(OUTPUT_DIR, "results.txt")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
cat("", file = RESULTS_LOG_PATH)

log_line <- function(...) {
  text <- paste0(...)
  cat(text, "\n", sep = "")
  cat(text, "\n", file = RESULTS_LOG_PATH, sep = "", append = TRUE)
}

log_section <- function(title = NULL) {
  log_line(strrep("=", 78))

  if (!is.null(title)) {
    log_line(title)
    log_line(strrep("=", 78))
  }
}

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
arm_data <- read.csv(
  file.path(ANALYSIS_DIR, "arm_level_data.csv"),
  stringsAsFactors = FALSE
)

STRATEGY_LEVELS <- c(
  "No screening / usual care",
  "PSA + DRE",
  "PSA",
  "PSA -> MRI",
  "Biomarker -> MRI",
  "MRI-first"
)
arm_data$strategy <- factor(arm_data$strategy, levels = STRATEGY_LEVELS)
stopifnot(!any(is.na(arm_data$strategy)))

# Internal consistency: GG1 cancers + clinically significant cancers cannot
# exceed the total number of prostate cancers in any arm.
stopifnot(all(arm_data$gg1 + arm_data$cs <= arm_data$pca))

# Wide (trial-level) frame for the randomized comparisons.
# RRs always use ITT (randomized / invited) denominators - the randomized
# comparison must not be conditioned on being screened.
intervention_arms <- arm_data[arm_data$arm == "I", ]
comparator_arms <- arm_data[arm_data$arm == "C", ]

intervention_arms <- intervention_arms[order(intervention_arms$trial), ]
comparator_arms <- comparator_arms[order(comparator_arms$trial), ]
stopifnot(identical(intervention_arms$trial, comparator_arms$trial))

trial_data <- data.frame(
  trial = intervention_arms$trial,
  year = intervention_arms$year,
  family = intervention_arms$family,
  label = paste0(intervention_arms$trial, " ", intervention_arms$year),
  intervention_n = intervention_arms$n_itt,
  comparator_n = comparator_arms$n_itt,
  intervention_gg1 = intervention_arms$gg1,
  comparator_gg1 = comparator_arms$gg1,
  intervention_cs = intervention_arms$cs,
  comparator_cs = comparator_arms$cs,
  intervention_pca = intervention_arms$pca,
  comparator_pca = comparator_arms$pca,
  intervention_biopsies = intervention_arms$bx,
  comparator_biopsies = comparator_arms$bx,
  stringsAsFactors = FALSE
)
trial_data$family_label <- ifelse(
  trial_data$family == "A",
  "A. Screening invitation vs no screening / usual care",
  "B. Diagnostic / biopsy strategy (both arms screened)"
)
trial_order <- order(trial_data$family, trial_data$year, trial_data$trial)
trial_data <- trial_data[trial_order, ]

# ---------------------------------------------------------------------------
# Descriptives
# ---------------------------------------------------------------------------
log_section("TABLE 1. Extracted data, internal consistency check")
descriptive_columns <- c(
  "trial",
  "arm",
  "strategy",
  "n_itt",
  "n_screened",
  "gg1",
  "cs",
  "pca",
  "bx"
)
descriptive_table <- arm_data[, descriptive_columns]
descriptive_table$gg1_pct_screened <- round(
  100 * descriptive_table$gg1 / descriptive_table$n_screened,
  3
)
descriptive_table$cs_pct_screened <- round(
  100 * descriptive_table$cs / descriptive_table$n_screened,
  3
)
descriptive_table$gg1_of_pca_pct <- round(
  100 * descriptive_table$gg1 / descriptive_table$pca,
  1
)
descriptive_table$bx_pct_screened <- round(
  100 * descriptive_table$bx / descriptive_table$n_screened,
  2
)

print(descriptive_table, row.names = FALSE)
capture.output(
  print(descriptive_table, row.names = FALSE),
  file = RESULTS_LOG_PATH,
  append = TRUE
)
write.csv(
  descriptive_table,
  file.path(OUTPUT_DIR, "table1_descriptives.csv"),
  row.names = FALSE
)
