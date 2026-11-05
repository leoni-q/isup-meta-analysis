###############################################################################
# PART 1. COMPARATIVE ANALYSES - RISK RATIOS
###############################################################################
log_section(
  "PART 1. COMPARATIVE META-ANALYSES (RISK RATIOS)"
)
log_line("Denominators are randomized / invited N for detection outcomes.")
log_line("Grade-mix outcomes use all detected cancers as denominator.")
log_line("Biopsy-yield outcomes use outcome-specific biopsy denominators from")
log_line("the supplied CSV and require data for both randomized arms.")
log_line("Biopsy burden (1.7) is the count of biopsies per randomized/invited")
log_line("participant, and is restricted to trials reporting both arms.")
log_line("Invitation and diagnostic-strategy contrasts are fitted separately;")
log_line("no pooled estimate is calculated across incompatible contrast families.")
log_line("PROBASE is reported separately as a timing-of-screening contrast.")
log_line("CAP remains in invitation-family models; its precision is not adjusted")
log_line("for practice-level clustering. Excluding-CAP models are in Part 3.")

RR_OUTCOMES <- list(
  gg1_detection = list(
    section = "1.1",
    analysis_name = "GG1 detection per randomized/invited participant",
    intervention_event_column = "intervention_gg1",
    intervention_total_column = "intervention_n",
    comparator_event_column = "comparator_gg1",
    comparator_total_column = "comparator_n",
    file_stem = "fig1_rr_gg1_detection",
    left_label = "Lower in intervention",
    right_label = "Higher in intervention",
    x_ticks = c(0.2, 0.5, 1, 2, 5, 10)
  ),
  gg1_share = list(
    section = "1.2",
    analysis_name = "GG1 among all detected prostate cancers",
    intervention_event_column = "intervention_gg1",
    intervention_total_column = "intervention_pca",
    comparator_event_column = "comparator_gg1",
    comparator_total_column = "comparator_pca",
    file_stem = "fig2_rr_gg1_share_of_pca",
    left_label = "Lower GG1 share",
    right_label = "Higher GG1 share",
    x_ticks = c(0.2, 0.35, 0.5, 0.7, 1, 1.5, 2)
  ),
  gg2plus_detection = list(
    section = "1.3",
    analysis_name = "GG>=2 detection per randomized/invited participant",
    intervention_event_column = "intervention_cs",
    intervention_total_column = "intervention_n",
    comparator_event_column = "comparator_cs",
    comparator_total_column = "comparator_n",
    file_stem = "fig3_rr_gg2plus_detection",
    left_label = "Lower in intervention",
    right_label = "Higher in intervention",
    x_ticks = c(0.2, 0.5, 1, 2, 5, 10)
  ),
  gg2plus_share = list(
    section = "1.4",
    analysis_name = "GG>=2 among all detected prostate cancers",
    intervention_event_column = "intervention_cs",
    intervention_total_column = "intervention_pca",
    comparator_event_column = "comparator_cs",
    comparator_total_column = "comparator_pca",
    file_stem = "fig4_rr_gg2plus_share_of_pca",
    left_label = "Lower GG>=2 share",
    right_label = "Higher GG>=2 share",
    x_ticks = c(0.2, 0.35, 0.5, 0.7, 1, 1.5, 2)
  ),
  gg1_per_biopsy = list(
    section = "1.5",
    analysis_name = "GG1 among participants biopsied",
    intervention_event_column = "intervention_gg1",
    intervention_total_column = "intervention_gg1_bx_n",
    comparator_event_column = "comparator_gg1",
    comparator_total_column = "comparator_gg1_bx_n",
    exclude_trials = "STHLM3-MRI",
    file_stem = "fig5_rr_gg1_per_biopsy",
    left_label = "Lower GG1 biopsy yield",
    right_label = "Higher GG1 biopsy yield",
    x_ticks = c(0.2, 0.35, 0.5, 0.7, 1, 1.5, 2)
  ),
  gg2plus_per_biopsy = list(
    section = "1.6",
    analysis_name = "GG>=2 among participants biopsied",
    intervention_event_column = "intervention_cs",
    intervention_total_column = "intervention_cs_bx_n",
    comparator_event_column = "comparator_cs",
    comparator_total_column = "comparator_cs_bx_n",
    exclude_trials = "STHLM3-MRI",
    file_stem = "fig6_rr_gg2plus_per_biopsy",
    left_label = "Lower GG>=2 biopsy yield",
    right_label = "Higher GG>=2 biopsy yield",
    x_ticks = c(0.2, 0.35, 0.5, 0.7, 1, 1.5, 2)
  ),
  # Protocol s.6 secondary endpoint: number of biopsies performed. This is the
  # randomized comparison of biopsy burden, distinct from the yield outcomes
  # above; only trials reporting biopsies in both arms contribute.
  biopsies = list(
    section = "1.7",
    analysis_name = "Biopsies performed per randomized/invited participant",
    intervention_event_column = "intervention_biopsies",
    intervention_total_column = "intervention_n",
    comparator_event_column = "comparator_biopsies",
    comparator_total_column = "comparator_n",
    file_stem = "fig7_rr_biopsies_per_participant",
    left_label = "Fewer biopsies",
    right_label = "More biopsies",
    x_ticks = c(0.25, 0.5, 1, 2, 4)
  )
)

risk_ratio_models <- list()

for (outcome_name in names(RR_OUTCOMES)) {
  outcome <- RR_OUTCOMES[[outcome_name]]

  for (cohort_name in names(RR_COHORTS)) {
    cohort <- RR_COHORTS[[cohort_name]]
    excluded_trials <- if (is.null(outcome$exclude_trials)) {
      character(0)
    } else {
      outcome$exclude_trials
    }
    is_excluded_trial <- cohort$data$trial %in% excluded_trials
    cohort_data <- cohort$data[!is_excluded_trial, ]
    if (any(is_excluded_trial)) {
      log_line(
        sprintf(
          "Excluded from %s (outcome-specific rule): %s",
          outcome$analysis_name,
          paste(unique(cohort$data$trial[is_excluded_trial]), collapse = ", ")
        )
      )
    }
    complete_cohort_data <- complete_rr_data(
      cohort_data,
      outcome$intervention_event_column,
      outcome$intervention_total_column,
      outcome$comparator_event_column,
      outcome$comparator_total_column
    )
    if (nrow(complete_cohort_data) == 0) {
      log_line(
        sprintf(
          "Skipping %s - %s: no trial reports both arms.",
          outcome$analysis_name,
          cohort$label
        )
      )
      next
    }

    model_title <- paste(outcome$analysis_name, cohort$label, sep = " - ")
    model <- fit_risk_ratio(
      source_data = cohort_data,
      intervention_event_column = outcome$intervention_event_column,
      intervention_total_column = outcome$intervention_total_column,
      comparator_event_column = outcome$comparator_event_column,
      comparator_total_column = outcome$comparator_total_column,
      title = model_title,
      include_subgroups = FALSE
    )
    model_key <- paste(outcome_name, cohort_name, sep = "__")
    risk_ratio_models[[model_key]] <- model

    report_risk_ratio(
      model,
      paste(
        outcome$section,
        toupper(outcome$analysis_name),
        "-",
        cohort$label
      )
    )
    collect_risk_ratio_results(
      model,
      outcome$analysis_name,
      cohort$label
    )
    report_leave_one_out_risk_ratio(
      source_data = cohort_data,
      intervention_event_column = outcome$intervention_event_column,
      intervention_total_column = outcome$intervention_total_column,
      comparator_event_column = outcome$comparator_event_column,
      comparator_total_column = outcome$comparator_total_column,
      analysis_name = outcome$analysis_name,
      group_name = cohort$label
    )
    save_risk_ratio_forest_plot(
      model,
      paste0(outcome$file_stem, "_", cohort_name),
      outcome$left_label,
      outcome$right_label,
      x_ticks = outcome$x_ticks
    )
  }

  complete_timing_data <- complete_rr_data(
    timing_trial_data,
    outcome$intervention_event_column,
    outcome$intervention_total_column,
    outcome$comparator_event_column,
    outcome$comparator_total_column
  )
  if (nrow(complete_timing_data) == 0) {
    log_line(
      sprintf(
        "Skipping %s - %s: both-arm data unavailable.",
        outcome$analysis_name,
        FAMILY_LABELS[["T"]]
      )
    )
    next
  }

  timing_model <- fit_risk_ratio(
    source_data = timing_trial_data,
    intervention_event_column = outcome$intervention_event_column,
    intervention_total_column = outcome$intervention_total_column,
    comparator_event_column = outcome$comparator_event_column,
    comparator_total_column = outcome$comparator_total_column,
    title = paste(outcome$analysis_name, FAMILY_LABELS[["T"]], sep = " - "),
    include_subgroups = FALSE
  )
  risk_ratio_models[[paste(outcome_name, "timing", sep = "__")]] <- timing_model
  report_risk_ratio(
    timing_model,
    paste(
      outcome$section,
      toupper(outcome$analysis_name),
      "- PROBASE STANDALONE TIMING CONTRAST"
    )
  )
  collect_risk_ratio_results(
    timing_model,
    outcome$analysis_name,
    FAMILY_LABELS[["T"]]
  )
  save_risk_ratio_forest_plot(
    timing_model,
    paste0(outcome$file_stem, "_timing"),
    outcome$left_label,
    outcome$right_label,
    x_ticks = outcome$x_ticks
  )
}
