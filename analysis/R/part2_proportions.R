###############################################################################
# PART 2. STUDY-LEVEL POOLED PROPORTIONS BY SCREENING STRATEGY
###############################################################################
log_section("PART 2. POOLED PROPORTIONS BY SCREENING STRATEGY (GLMM)")
log_line("Each trial contributes at most one independent unit to a strategy.")
log_line("Same-strategy arms within a trial are combined before model fitting.")
log_line("PROBASE contributes only its immediate-PSA arm.")
log_line("GG1 and GG>=2 detection rates use the full randomized/invited arm.")
log_line("Grade mix uses all detected prostate cancers as denominator.")
log_line("Biopsy yield uses outcome-specific biopsy denominators from the CSV.")
log_line("Strategies are fitted separately; no all-strategy estimate or test is run.")
log_line("")
log_line("Study-strategy units:")

for (strategy in levels(droplevels(proportion_data$strategy))) {
  uses_strategy <- proportion_data$strategy == strategy
  study_labels <- paste0(
    proportion_data$trial[uses_strategy],
    " (",
    proportion_data$source_arms[uses_strategy],
    ")"
  )
  log_line(
    sprintf(
      "  %-28s %s",
      strategy,
      paste(study_labels, collapse = ", ")
    )
  )
}

PROPORTION_OUTCOMES <- list(
  gg1_detection = list(
    section = "2.1",
    analysis_name = "GG1 detection per randomized/invited participant",
    event_column = "gg1",
    total_column = "n_itt",
    file_stem = "fig7_prop_gg1_per_randomized",
    axis_label = "ISUP GG1 detected per 100 randomized/invited (%)"
  ),
  gg1_share = list(
    section = "2.2",
    analysis_name = "GG1 among all detected prostate cancers",
    event_column = "gg1",
    total_column = "pca",
    file_stem = "fig8_prop_gg1_share_of_pca",
    axis_label = "ISUP GG1 as % of all detected prostate cancers"
  ),
  gg2plus_detection = list(
    section = "2.3",
    analysis_name = "GG>=2 detection per randomized/invited participant",
    event_column = "cs",
    total_column = "n_itt",
    file_stem = "fig9_prop_gg2plus_per_randomized",
    axis_label = "ISUP GG>=2 detected per 100 randomized/invited (%)"
  ),
  gg2plus_share = list(
    section = "2.4",
    analysis_name = "GG>=2 among all detected prostate cancers",
    event_column = "cs",
    total_column = "pca",
    file_stem = "fig10_prop_gg2plus_share_of_pca",
    axis_label = "ISUP GG>=2 as % of all detected prostate cancers"
  ),
  gg1_per_biopsy = list(
    section = "2.5",
    analysis_name = "GG1 among participants biopsied",
    event_column = "gg1",
    total_column = "gg1_bx_n",
    file_stem = "fig11_prop_gg1_per_biopsy",
    axis_label = "ISUP GG1 per 100 participants biopsied (%)"
  ),
  gg2plus_per_biopsy = list(
    section = "2.6",
    analysis_name = "GG>=2 among participants biopsied",
    event_column = "cs",
    total_column = "cs_bx_n",
    file_stem = "fig12_prop_gg2plus_per_biopsy",
    axis_label = "ISUP GG>=2 per 100 participants biopsied (%)"
  ),
  biopsies = list(
    section = "2.7",
    analysis_name = "Biopsies per screened/assessed participant",
    event_column = "bx",
    total_column = "n_screened",
    file_stem = "fig13_prop_biopsies_per_screened",
    axis_label = "Biopsies per 100 screened/assessed (%)"
  )
)

proportion_models <- list()

strategy_slug <- function(strategy) {
  slug <- gsub("[^a-z0-9]+", "_", tolower(strategy))
  gsub("^_|_$", "", slug)
}

for (outcome_name in names(PROPORTION_OUTCOMES)) {
  outcome <- PROPORTION_OUTCOMES[[outcome_name]]
  is_complete <- (
    !is.na(proportion_data[[outcome$event_column]]) &
      !is.na(proportion_data[[outcome$total_column]])
  )
  outcome_data <- proportion_data[is_complete, ]
  included_strategies <- levels(droplevels(outcome_data$strategy))

  for (strategy in included_strategies) {
    strategy_data <- outcome_data[outcome_data$strategy == strategy, ]
    model <- fit_proportion(
      source_data = strategy_data,
      event_column = outcome$event_column,
      total_column = outcome$total_column,
      title = paste(outcome$analysis_name, strategy, sep = " - "),
      include_subgroups = FALSE
    )
    model_key <- paste(outcome_name, strategy_slug(strategy), sep = "__")
    proportion_models[[model_key]] <- model

    report_proportion(
      model,
      paste(outcome$section, toupper(outcome$analysis_name), "-", strategy),
      strategy
    )
    collect_proportion_results(model, outcome$analysis_name, strategy)
    report_leave_one_out_proportion(
      source_data = strategy_data,
      event_column = outcome$event_column,
      total_column = outcome$total_column,
      analysis_name = outcome$analysis_name,
      group_name = strategy
    )
    save_proportion_forest_plot(
      model,
      paste(outcome$file_stem, strategy_slug(strategy), sep = "_"),
      outcome$axis_label,
      summary_label = paste("Pooled", strategy, "(%)")
    )
  }
}
