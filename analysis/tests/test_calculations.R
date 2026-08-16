# Unit and reproducibility tests for the meta-analysis calculation helpers.
#
# Run from anywhere:
#   Rscript analysis/tests/test_calculations.R
#
# No testthat dependency: each check is registered through expect() and the
# script exits non-zero if any check fails.
#
# Three layers:
#   1. Unit tests of the pure helpers (fit_risk_ratio, complete_rr_data,
#      extract_proportion_result, format_heterogeneity, fit_proportion)
#      against hand-computed values and independent metafor refits.
#   2. Integration test of the data preparation in setup_and_data.R,
#      run in a temporary sandbox so no real outputs are touched.
#   3. Reproducibility test: headline models are refitted from the CSV and
#      compared against the published table2_pooled_estimates.csv.

suppressPackageStartupMessages({
  library(meta)
  library(metafor)
})

# --- tiny test harness -------------------------------------------------------
test_failures <- character(0)
test_count <- 0

expect <- function(description, condition) {
  test_count <<- test_count + 1

  ok <- isTRUE(condition)
  if (!ok) {
    test_failures <<- c(test_failures, description)
  }
  cat(sprintf("  [%s] %s\n", if (ok) "PASS" else "FAIL", description))
  invisible(ok)
}

expect_close <- function(description, actual, expected, tolerance = 1e-8) {
  expect(
    sprintf("%s (actual %.6g, expected %.6g)", description, actual, expected),
    is.finite(actual) && is.finite(expected) &&
      abs(actual - expected) <= tolerance
  )
}

command_arguments <- commandArgs(trailingOnly = FALSE)
script_argument <- command_arguments[grepl("^--file=", command_arguments)]
TEST_DIR <- dirname(normalizePath(sub("^--file=", "", script_argument[1])))
REAL_ANALYSIS_DIR <- dirname(TEST_DIR)
MODULE_DIR <- file.path(REAL_ANALYSIS_DIR, "R")

# Helper modules use log_line/log_section only inside report functions, which
# these tests do not call, but stubs keep any accidental call harmless.
log_line <- function(...) invisible(NULL)
log_section <- function(...) invisible(NULL)

source(file.path(MODULE_DIR, "summary_helpers.R"))
source(file.path(MODULE_DIR, "risk_ratio_helpers.R"))
source(file.path(MODULE_DIR, "proportion_helpers.R"))

###############################################################################
cat("1. complete_rr_data\n")
###############################################################################
mixed_data <- data.frame(
  trial = c("T1", "T2", "T3", "T4"),
  label = c("T1", "T2", "T3", "T4"),
  e_events = c(10, NA, 12, 13),
  e_total = c(100, 100, NA, 100),
  c_events = c(5, 6, 7, 8),
  c_total = c(100, 100, 100, NA)
)
kept <- complete_rr_data(mixed_data, "e_events", "e_total", "c_events", "c_total")
expect(
  "keeps only rows with all four cells non-missing",
  identical(kept$trial, "T1")
)
expect(
  "returns a zero-row frame when nothing is complete",
  nrow(complete_rr_data(
    mixed_data[2:4, ], "e_events", "e_total", "c_events", "c_total"
  )) == 0
)

###############################################################################
cat("2. fit_risk_ratio: single study, no zero cells\n")
###############################################################################
single_study <- data.frame(
  trial = "T1", label = "T1",
  e_events = 30, e_total = 200, c_events = 10, c_total = 190
)
single_model <- fit_risk_ratio(
  single_study, "e_events", "e_total", "c_events", "c_total",
  title = "test", include_subgroups = FALSE
)
expect_close(
  "point estimate equals the raw risk ratio (a/n1)/(c/n2)",
  exp(single_model$TE.random),
  (30 / 200) / (10 / 190)
)
expect_close(
  "log-RR variance equals 1/a - 1/n1 + 1/c - 1/n2",
  single_model$seTE[1]^2,
  1 / 30 - 1 / 200 + 1 / 10 - 1 / 190
)

###############################################################################
cat("3. fit_risk_ratio: DL tau and Hartung-Knapp CI match metafor\n")
###############################################################################
two_studies <- data.frame(
  trial = c("T1", "T2"), label = c("T1", "T2"),
  e_events = c(30, 45), e_total = c(200, 400),
  c_events = c(10, 40), c_total = c(190, 410)
)
meta_model <- fit_risk_ratio(
  two_studies, "e_events", "e_total", "c_events", "c_total",
  title = "test", include_subgroups = FALSE
)
log_rr <- log((two_studies$e_events / two_studies$e_total) /
  (two_studies$c_events / two_studies$c_total))
log_rr_var <- 1 / two_studies$e_events - 1 / two_studies$e_total +
  1 / two_studies$c_events - 1 / two_studies$c_total
metafor_model <- rma(
  yi = log_rr, vi = log_rr_var, method = "DL", test = "knha"
)
expect_close(
  "pooled log RR matches independent metafor DL fit",
  meta_model$TE.random, as.numeric(metafor_model$beta)
)
expect_close(
  "DL tau^2 matches metafor",
  meta_model$tau2, metafor_model$tau2
)
expect_close(
  "Hartung-Knapp lower CI matches metafor test='knha'",
  meta_model$lower.random, metafor_model$ci.lb
)
expect_close(
  "Hartung-Knapp upper CI matches metafor test='knha'",
  meta_model$upper.random, metafor_model$ci.ub
)
expect_close(
  "Cochran's Q matches metafor",
  meta_model$Q, metafor_model$QE
)

###############################################################################
cat("4. fit_risk_ratio: zero-cell continuity handling (protocol: incr 0.5)\n")
###############################################################################
zero_cell_study <- data.frame(
  trial = "T1", label = "T1",
  e_events = 33, e_total = 23341, c_events = 0, c_total = 23301
)
zero_cell_model <- fit_risk_ratio(
  zero_cell_study, "e_events", "e_total", "c_events", "c_total",
  title = "test", include_subgroups = FALSE
)
expect(
  "zero-cell RR is finite after 0.5 continuity correction",
  is.finite(exp(zero_cell_model$TE.random))
)
# meta with RR.Cochrane = TRUE, incr = 0.5, method.incr = "only0" computes
# ((a+0.5)/(n1+1)) / ((c+0.5)/(n2+1)) for the zero-cell study.
expect_close(
  "zero-cell RR matches the Cochrane 0.5-increment formula",
  exp(zero_cell_model$TE.random),
  ((33 + 0.5) / (23341 + 1)) / ((0 + 0.5) / (23301 + 1)),
  tolerance = 1e-6
)
no_zero_model <- fit_risk_ratio(
  single_study, "e_events", "e_total", "c_events", "c_total",
  title = "test", include_subgroups = FALSE
)
expect_close(
  "continuity correction is NOT applied when no cell is zero (only0)",
  exp(no_zero_model$TE.random),
  (30 / 200) / (10 / 190)
)

###############################################################################
cat("5. extract_proportion_result\n")
###############################################################################
one_study_prop <- data.frame(
  trial = "T1", year = 2020, arm = "I",
  events = 30, total = 400,
  strategy = factor("PSA")
)
one_prop_model <- fit_proportion(
  one_study_prop, "events", "total", title = "test"
)
exact_ci <- binom.test(30, 400)$conf.int
one_prop_result <- extract_proportion_result(one_prop_model)
expect_close(
  "k=1 estimate is the raw percentage",
  one_prop_result$estimate, 100 * 30 / 400
)
expect_close(
  "k=1 lower bound is the exact (Clopper-Pearson) limit",
  one_prop_result$lower, 100 * exact_ci[1]
)
expect_close(
  "k=1 upper bound is the exact (Clopper-Pearson) limit",
  one_prop_result$upper, 100 * exact_ci[2]
)
expect(
  "k=1 heterogeneity fields are NA",
  is.na(one_prop_result$Q) && is.na(one_prop_result$I2)
)

multi_study_prop <- data.frame(
  trial = c("T1", "T2", "T3"), year = c(2020, 2021, 2022),
  arm = c("I", "I", "I"),
  events = c(30, 55, 12), total = c(400, 600, 300),
  strategy = factor(c("PSA", "PSA", "PSA"))
)
multi_prop_model <- fit_proportion(
  multi_study_prop, "events", "total", title = "test"
)
multi_prop_result <- extract_proportion_result(multi_prop_model)
expect_close(
  "k>1 estimate is the back-transformed (plogis) GLMM logit",
  multi_prop_result$estimate, 100 * plogis(multi_prop_model$TE.random)
)
expect_close(
  "k>1 CI bounds are back-transformed with plogis",
  multi_prop_result$lower, 100 * plogis(multi_prop_model$lower.random)
)
observed_percent <- 100 * multi_study_prop$events / multi_study_prop$total
expect(
  "pooled GLMM proportion lies within the observed study range",
  multi_prop_result$estimate >= min(observed_percent) &&
    multi_prop_result$estimate <= max(observed_percent)
)
expect(
  "GLMM model is fitted with ML tau and PLOGIT transform",
  multi_prop_model$method == "GLMM" &&
    multi_prop_model$method.tau == "ML" &&
    multi_prop_model$sm == "PLOGIT"
)

###############################################################################
cat("6. format_heterogeneity\n")
###############################################################################
heterogeneity_text <- format_heterogeneity(meta_model)
expect(
  "reports Cochran's Q, I2 and tau2",
  grepl("Cochran's Q", heterogeneity_text) &&
    grepl("I2", heterogeneity_text) &&
    grepl("tau2", heterogeneity_text)
)
expect(
  "flags significance only when p(Q) < 0.05",
  grepl("significant heterogeneity", heterogeneity_text) ==
    (meta_model$pval.Q[1] < 0.05)
)

###############################################################################
cat("7. Protocol settings on the fitted RR model object\n")
###############################################################################
expect(
  "RR model uses inverse variance, DL tau, Hartung-Knapp, incr 0.5/only0",
  meta_model$method == "Inverse" &&
    meta_model$method.tau == "DL" &&
    identical(unname(meta_model$method.random.ci), "HK") &&
    meta_model$incr == 0.5 &&
    meta_model$method.incr == "only0" &&
    meta_model$sm == "RR"
)

###############################################################################
cat("8. Data preparation (setup_and_data.R) in a sandbox\n")
###############################################################################
# setup_and_data.R deletes figures and writes tables in ANALYSIS_DIR/output,
# so run it against a temporary copy of the CSV, never the real output dir.
sandbox_dir <- file.path(tempdir(), "meta_analysis_test_sandbox")
dir.create(sandbox_dir, showWarnings = FALSE, recursive = TRUE)
invisible(file.copy(
  file.path(REAL_ANALYSIS_DIR, "arm_level_data.csv"),
  file.path(sandbox_dir, "arm_level_data.csv"),
  overwrite = TRUE
))
ANALYSIS_DIR <- sandbox_dir
setup_output <- capture.output(
  source(file.path(MODULE_DIR, "setup_and_data.R"))
)

erspc_row <- trial_data[trial_data$trial == "ERSPC", ]
expect(
  "intervention/comparator arms are paired without swapping (ERSPC)",
  erspc_row$intervention_gg1 == 4325 && erspc_row$comparator_gg1 == 2360 &&
    erspc_row$intervention_n == 72890 && erspc_row$comparator_n == 89353
)
expect(
  "RR cohorts split by family: A=4, B=4, T=PROBASE only",
  nrow(invitation_trial_data) == 4 &&
    nrow(diagnostic_trial_data) == 4 &&
    identical(timing_trial_data$trial, "PROBASE")
)

goteborg_unit <- proportion_data[proportion_data$trial == "GOTEBORG-2", ]
expect(
  "GOTEBORG-2 same-strategy arms are summed (events and denominators)",
  nrow(goteborg_unit) == 1 &&
    goteborg_unit$gg1 == 68 + 159 &&
    goteborg_unit$cs == 117 + 139 &&
    goteborg_unit$pca == 185 + 298 &&
    goteborg_unit$n_itt == 6575 + 6578 &&
    goteborg_unit$bx == 254 + 588
)
probase_units <- proportion_data[proportion_data$trial == "PROBASE", ]
expect(
  "PROBASE contributes only its immediate-PSA (I) arm to proportions",
  nrow(probase_units) == 1 && probase_units$source_arms == "I"
)
expect(
  "sum_or_na: all-NA collapses to NA (ERSPC biopsy count stays missing)",
  is.na(proportion_data$bx[
    proportion_data$trial == "ERSPC" &
      proportion_data$strategy == "No screening / usual care"
  ])
)

###############################################################################
cat("9. Reproducibility: refits match published table2 values\n")
###############################################################################
published <- read.csv(
  file.path(REAL_ANALYSIS_DIR, "output", "table2_pooled_estimates.csv"),
  stringsAsFactors = FALSE
)
published_row <- function(analysis_name, group_prefix) {
  rows <- published[
    published$analysis == analysis_name &
      startsWith(published$group, group_prefix),
  ]
  stopifnot(nrow(rows) == 1)
  rows
}

refit_gg1_invitation <- fit_risk_ratio(
  invitation_trial_data,
  "intervention_gg1", "intervention_n", "comparator_gg1", "comparator_n",
  title = "test", include_subgroups = FALSE
)
gg1_row <- published_row(
  "GG1 detection per randomized/invited participant", "A."
)
expect_close(
  "GG1 detection RR (invitation family) reproduces table2",
  exp(refit_gg1_invitation$TE.random), gg1_row$estimate,
  tolerance = 1e-3
)
expect_close(
  "GG1 detection RR lower CI reproduces table2",
  exp(refit_gg1_invitation$lower.random), gg1_row$ci_lower,
  tolerance = 1e-3
)

refit_probase_cs <- fit_risk_ratio(
  timing_trial_data,
  "intervention_cs", "intervention_n", "comparator_cs", "comparator_n",
  title = "test", include_subgroups = FALSE
)
probase_row <- published_row(
  "GG>=2 detection per randomized/invited participant", "T."
)
expect_close(
  "PROBASE GG>=2 detection RR (zero-cell arm) reproduces table2",
  exp(refit_probase_cs$TE.random), probase_row$estimate,
  tolerance = 1e-3
)

refit_gg2_diag_bx <- fit_risk_ratio(
  diagnostic_trial_data,
  "intervention_cs", "intervention_cs_bx_n",
  "comparator_cs", "comparator_cs_bx_n",
  title = "test", include_subgroups = FALSE
)
gg2_bx_row <- published_row("GG>=2 among participants biopsied", "B.")
expect_close(
  "GG>=2 biopsy-yield RR (diagnostic family) reproduces table2",
  exp(refit_gg2_diag_bx$TE.random), gg2_bx_row$estimate,
  tolerance = 1e-3
)

unlink(sandbox_dir, recursive = TRUE)

###############################################################################
cat(sprintf(
  "\n%d checks run, %d failed.\n", test_count, length(test_failures)
))
if (length(test_failures) > 0) {
  cat("Failed checks:\n")
  cat(paste0("  - ", test_failures, "\n"), sep = "")
  quit(status = 1)
}
cat("All calculation tests passed.\n")
