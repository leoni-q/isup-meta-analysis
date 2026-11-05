# Package setup, output paths, logging, and data preparation.

suppressPackageStartupMessages({
  library(meta)
  library(metafor)
})

settings.meta(digits = 3, digits.pval = 4)

OUTPUT_DIR <- file.path(ANALYSIS_DIR, "output")
RESULTS_LOG_PATH <- file.path(OUTPUT_DIR, "results.txt")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Clearing previous outputs is deferred until the data checks below have passed
# (see clear_previous_outputs()), so an aborted run leaves the last good
# figures and results.txt intact instead of an empty output directory.
clear_previous_outputs <- function() {
  obsolete_plot_paths <- list.files(
    OUTPUT_DIR,
    pattern = "^fig.*\\.(pdf|png)$",
    full.names = TRUE
  )
  if (length(obsolete_plot_paths) > 0) {
    unlink(obsolete_plot_paths)
  }
  cat("", file = RESULTS_LOG_PATH)
}

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

# ProScreen is excluded from every analysis. Its rows stay in the CSV so the
# extraction remains complete, and are dropped here before any cohort,
# study-strategy unit or descriptive table is built.
EXCLUDED_TRIALS <- c("ProScreen")
stopifnot(all(EXCLUDED_TRIALS %in% arm_data$trial))
arm_data <- arm_data[!arm_data$trial %in% EXCLUDED_TRIALS, ]

STRATEGY_LEVELS <- c(
  "No screening / usual care",
  "Delayed screening / DRE subset",
  "PSA + DRE",
  "PSA",
  "PSA -> MRI",
  "Biomarker -> MRI",
  "MRI-first"
)
arm_data$strategy <- factor(arm_data$strategy, levels = STRATEGY_LEVELS)
stopifnot(!any(is.na(arm_data$strategy)))

# Internal consistency: GG1 cancers + clinically significant cancers cannot
# exceed the total number of prostate cancers in any arm. Biopsy-denominator
# fields are kept separately because the supplied CSV marks some arm/outcome
# combinations unavailable even when a biopsy count exists (notably mature CAP
# cancer counts versus first-round biopsies).
stopifnot(
  all(arm_data$gg1 + arm_data$cs <= arm_data$pca),
  all(arm_data$gg1[!is.na(arm_data$gg1_bx_n)] <=
    arm_data$gg1_bx_n[!is.na(arm_data$gg1_bx_n)]),
  all(arm_data$cs[!is.na(arm_data$cs_bx_n)] <=
    arm_data$cs_bx_n[!is.na(arm_data$cs_bx_n)])
)

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
  intervention_n_screened = intervention_arms$n_screened,
  comparator_n_screened = comparator_arms$n_screened,
  intervention_gg1 = intervention_arms$gg1,
  comparator_gg1 = comparator_arms$gg1,
  intervention_cs = intervention_arms$cs,
  comparator_cs = comparator_arms$cs,
  intervention_pca = intervention_arms$pca,
  comparator_pca = comparator_arms$pca,
  intervention_biopsies = intervention_arms$bx,
  comparator_biopsies = comparator_arms$bx,
  intervention_gg1_bx_n = intervention_arms$gg1_bx_n,
  comparator_gg1_bx_n = comparator_arms$gg1_bx_n,
  intervention_cs_bx_n = intervention_arms$cs_bx_n,
  comparator_cs_bx_n = comparator_arms$cs_bx_n,
  stringsAsFactors = FALSE
)
FAMILY_LABELS <- c(
  A = "A. Screening invitation vs no screening / usual care",
  B = "B. Diagnostic / biopsy strategy (both arms screened)",
  T = "T. Timing of screening initiation"
)
stopifnot(all(trial_data$family %in% names(FAMILY_LABELS)))
trial_data$family_label <- unname(FAMILY_LABELS[trial_data$family])
trial_order <- order(trial_data$family, trial_data$year, trial_data$trial)
trial_data <- trial_data[trial_order, ]

# Explicit analysis cohorts. PROBASE is a timing contrast and is never pooled
# with invitation or diagnostic-strategy trials. CAP remains in the primary
# invitation cohort; its practice-level clustering is documented and tested by
# an excluding-CAP sensitivity analysis.
invitation_trial_data <- trial_data[trial_data$family == "A", ]
diagnostic_trial_data <- trial_data[trial_data$family == "B", ]
timing_trial_data <- trial_data[trial_data$family == "T", ]

stopifnot(
  identical(
    sort(invitation_trial_data$trial),
    sort(c("CAP", "ERSPC", "PLCO"))
  ),
  identical(
    sort(diagnostic_trial_data$trial),
    sort(c("GOTEBORG-2", "MVP", "PROSA", "STHLM3-MRI"))
  ),
  identical(timing_trial_data$trial, "PROBASE")
)

RR_COHORTS <- list(
  invitation = list(
    label = FAMILY_LABELS[["A"]],
    data = invitation_trial_data
  ),
  diagnostic = list(
    label = FAMILY_LABELS[["B"]],
    data = diagnostic_trial_data
  )
)

# Proportion models use one independent contribution per trial and strategy.
# PROBASE contributes only its immediate-PSA arm because the delayed arm had
# not undergone PSA screening at this endpoint. If both randomized arms use
# the same strategy (GOTEBORG-2), counts and denominators are summed into one
# study-strategy contribution.
proportion_arm_data <- arm_data[
  arm_data$trial != "PROBASE" | arm_data$arm == "I",
]

sum_or_na <- function(values) {
  if (all(is.na(values))) {
    return(NA_real_)
  }

  sum(values, na.rm = TRUE)
}

proportion_group_key <- interaction(
  proportion_arm_data$trial,
  proportion_arm_data$strategy,
  drop = TRUE,
  lex.order = TRUE
)
proportion_groups <- split(proportion_arm_data, proportion_group_key)

proportion_data <- do.call(
  rbind,
  lapply(proportion_groups, function(group_data) {
    data.frame(
      trial = group_data$trial[1],
      year = group_data$year[1],
      strategy = as.character(group_data$strategy[1]),
      source_arms = paste(group_data$arm, collapse = "+"),
      arm_count = nrow(group_data),
      n_itt = sum_or_na(group_data$n_itt),
      n_screened = sum_or_na(group_data$n_screened),
      gg1 = sum_or_na(group_data$gg1),
      cs = sum_or_na(group_data$cs),
      pca = sum_or_na(group_data$pca),
      bx = sum_or_na(group_data$bx),
      gg1_bx_n = sum_or_na(group_data$gg1_bx_n),
      cs_bx_n = sum_or_na(group_data$cs_bx_n),
      stringsAsFactors = FALSE
    )
  })
)
row.names(proportion_data) <- NULL
proportion_data$strategy <- factor(
  proportion_data$strategy,
  levels = STRATEGY_LEVELS
)
proportion_data <- proportion_data[
  order(proportion_data$strategy, proportion_data$year, proportion_data$trial),
]

stopifnot(
  !any(duplicated(proportion_data[c("trial", "strategy")])),
  nrow(proportion_data[
    proportion_data$trial == "PROBASE" &
      proportion_data$strategy == "PSA",
  ]) == 1,
  proportion_data$source_arms[
    proportion_data$trial == "PROBASE"
  ] == "I",
  proportion_data$arm_count[
    proportion_data$trial == "GOTEBORG-2"
  ] == 2,
  all(proportion_data$gg1 + proportion_data$cs <= proportion_data$pca),
  all(proportion_data$gg1 <= proportion_data$n_itt),
  all(proportion_data$cs <= proportion_data$n_itt),
  all(proportion_data$gg1[!is.na(proportion_data$gg1_bx_n)] <=
    proportion_data$gg1_bx_n[!is.na(proportion_data$gg1_bx_n)]),
  all(proportion_data$cs[!is.na(proportion_data$cs_bx_n)] <=
    proportion_data$cs_bx_n[!is.na(proportion_data$cs_bx_n)])
)

# Every structural check on the extracted data has now passed, so replacing the
# previous run's outputs is safe.
clear_previous_outputs()

# ---------------------------------------------------------------------------
# Descriptives
# ---------------------------------------------------------------------------
log_section("TABLE 1. Extracted data, internal consistency check")
log_line(
  sprintf(
    "Excluded from all analyses: %s",
    paste(EXCLUDED_TRIALS, collapse = ", ")
  )
)
descriptive_columns <- c(
  "trial",
  "arm",
  "strategy",
  "n_itt",
  "n_screened",
  "n_screened_source",
  "gg1",
  "cs",
  "pca",
  "bx",
  "gg1_bx_n",
  "cs_bx_n",
  "flag"
)
descriptive_table <- arm_data[, descriptive_columns]
descriptive_table$gg1_pct_itt <- round(
  100 * descriptive_table$gg1 / descriptive_table$n_itt,
  3
)
descriptive_table$cs_pct_itt <- round(
  100 * descriptive_table$cs / descriptive_table$n_itt,
  3
)
descriptive_table$gg1_of_pca_pct <- round(
  100 * descriptive_table$gg1 / descriptive_table$pca,
  1
)
descriptive_table$cs_of_pca_pct <- round(
  100 * descriptive_table$cs / descriptive_table$pca,
  1
)
descriptive_table$gg1_of_bx_pct <- round(
  100 * descriptive_table$gg1 / descriptive_table$gg1_bx_n,
  1
)
descriptive_table$cs_of_bx_pct <- round(
  100 * descriptive_table$cs / descriptive_table$cs_bx_n,
  1
)
descriptive_table$bx_pct_screened <- round(
  100 * descriptive_table$bx / descriptive_table$n_screened,
  2
)

print(descriptive_table, row.names = FALSE)
descriptive_output <- capture.output(
  print(descriptive_table, row.names = FALSE)
)
descriptive_output <- sub("[[:space:]]+$", "", descriptive_output)
cat(
  paste0(descriptive_output, "\n"),
  file = RESULTS_LOG_PATH,
  append = TRUE,
  sep = ""
)
write.csv(
  descriptive_table,
  file.path(OUTPUT_DIR, "table1_descriptives.csv"),
  row.names = FALSE
)
