###############################################################################
# PART 3. SENSITIVITY ANALYSES
###############################################################################
log_section("PART 3. SENSITIVITY ANALYSES")

# 3.1 Proportions on ITT denominators
log_line("")
log_line("--- 3.1 Proportions on ITT (randomized/invited) denominators ---")

gg1_itt_proportion <- fit_proportion(
  source_data = arm_data,
  event_column = "gg1",
  total_column = "n_itt",
  title = "GG1 per man randomized/invited, by strategy"
)
report_proportion(
  gg1_itt_proportion,
  "3.1a GG1 detection rate per man randomized/invited, by strategy"
)
collect_proportion_results(
  gg1_itt_proportion,
  "SENS: GG1 rate per man randomized (proportion)"
)
save_proportion_forest_plot(
  gg1_itt_proportion,
  "figS1_prop_gg1_per_randomized",
  "ISUP GG1 detected per 100 men randomized/invited (%)"
)

cspca_itt_proportion <- fit_proportion(
  source_data = arm_data,
  event_column = "cs",
  total_column = "n_itt",
  title = "csPCa per man randomized/invited, by strategy"
)
report_proportion(
  cspca_itt_proportion,
  "3.1b csPCa detection rate per man randomized/invited, by strategy"
)
collect_proportion_results(
  cspca_itt_proportion,
  "SENS: csPCa rate per man randomized (proportion)"
)
save_proportion_forest_plot(
  cspca_itt_proportion,
  "figS2_prop_cspca_per_randomized",
  "ISUP >=2 detected per 100 men randomized/invited (%)"
)

# Each list item is a small object describing one outcome. The same outcome
# definitions are reused across the trial-exclusion and denominator analyses.
DETECTION_OUTCOMES <- list(
  gg1 = list(
    intervention_event_column = "intervention_gg1",
    comparator_event_column = "comparator_gg1",
    label = "GG1 detection"
  ),
  cspca = list(
    intervention_event_column = "intervention_cs",
    comparator_event_column = "comparator_cs",
    label = "csPCa detection"
  ),
  any_pca = list(
    intervention_event_column = "intervention_pca",
    comparator_event_column = "comparator_pca",
    label = "any PCa detection"
  )
)

# 3.2 Exclude PROBASE (delayed screening comparator; control csPCa = 0 inferred)
log_line("")
log_line(
  paste0(
    "--- 3.2 RRs excluding PROBASE (delayed-screening comparator; ",
    "control csPCa = 0 inferred) ---"
  )
)
trials_without_probase <- trial_data[trial_data$trial != "PROBASE", ]

for (outcome in DETECTION_OUTCOMES) {
  sensitivity_model <- fit_risk_ratio(
    source_data = trials_without_probase,
    intervention_event_column = outcome$intervention_event_column,
    intervention_total_column = "intervention_n",
    comparator_event_column = outcome$comparator_event_column,
    comparator_total_column = "comparator_n",
    title = paste0(outcome$label, " (PROBASE excluded)")
  )
  report_risk_ratio(
    sensitivity_model,
    paste0("3.2 ", outcome$label, " RR - PROBASE excluded")
  )
  collect_risk_ratio_results(
    sensitivity_model,
    paste0("SENS no PROBASE: ", outcome$label, " (RR)")
  )
}

gg1_share_without_probase <- fit_risk_ratio(
  source_data = trials_without_probase,
  intervention_event_column = "intervention_gg1",
  intervention_total_column = "intervention_pca",
  comparator_event_column = "comparator_gg1",
  comparator_total_column = "comparator_pca",
  title = "GG1 share of PCa (PROBASE excluded)"
)
report_risk_ratio(
  gg1_share_without_probase,
  "3.2 GG1 share of detected PCa RR - PROBASE excluded"
)
collect_risk_ratio_results(
  gg1_share_without_probase,
  "SENS no PROBASE: GG1 share of detected PCa (RR)"
)

# 3.3 Exclude ERSPC (Gleason counts back-calculated from percentages)
log_line("")
log_line(
  paste0(
    "--- 3.3 RRs excluding ERSPC ",
    "(Gleason counts back-calculated from reported percentages) ---"
  )
)
trials_without_erspc <- trial_data[trial_data$trial != "ERSPC", ]

for (outcome in DETECTION_OUTCOMES[c("gg1", "cspca")]) {
  sensitivity_model <- fit_risk_ratio(
    source_data = trials_without_erspc,
    intervention_event_column = outcome$intervention_event_column,
    intervention_total_column = "intervention_n",
    comparator_event_column = outcome$comparator_event_column,
    comparator_total_column = "comparator_n",
    title = paste0(outcome$label, " (ERSPC excluded)")
  )
  report_risk_ratio(
    sensitivity_model,
    paste0("3.3 ", outcome$label, " RR - ERSPC excluded")
  )
  collect_risk_ratio_results(
    sensitivity_model,
    paste0("SENS no ERSPC: ", outcome$label, " (RR)")
  )
}

# 3.4 MVP and PROSA on randomized rather than analysed denominators
log_line("")
log_line(
  "--- 3.4 RRs with MVP and PROSA on randomized (not analysed) denominators ---"
)
randomized_denominator_data <- trial_data
is_mvp_trial <- randomized_denominator_data$trial == "MVP"
is_prosa_trial <- randomized_denominator_data$trial == "PROSA"
randomized_denominator_data$intervention_n[is_mvp_trial] <- 259
randomized_denominator_data$comparator_n[is_mvp_trial] <- 266
randomized_denominator_data$intervention_n[is_prosa_trial] <- 374
randomized_denominator_data$comparator_n[is_prosa_trial] <- 385

for (outcome in DETECTION_OUTCOMES[c("gg1", "cspca")]) {
  sensitivity_model <- fit_risk_ratio(
    source_data = randomized_denominator_data,
    intervention_event_column = outcome$intervention_event_column,
    intervention_total_column = "intervention_n",
    comparator_event_column = outcome$comparator_event_column,
    comparator_total_column = "comparator_n",
    title = paste0(outcome$label, " (randomized denominators)")
  )
  report_risk_ratio(
    sensitivity_model,
    paste0(
      "3.4 ",
      outcome$label,
      " RR - MVP/PROSA randomized denominators"
    )
  )
  collect_risk_ratio_results(
    sensitivity_model,
    paste0("SENS randomized denom: ", outcome$label, " (RR)")
  )
}

# 3.5 Continuity-correction sensitivity for the single zero cell (PROBASE csPCa)
log_line("")
log_line(
  paste0(
    "--- 3.5 Continuity-correction sensitivity, csPCa RR ",
    "(single zero cell: PROBASE control) ---"
  )
)

for (continuity_correction in c(0.1, 0.5, 1)) {
  continuity_model <- suppressWarnings(
    metabin(
      event.e = trial_data$intervention_cs,
      n.e = trial_data$intervention_n,
      event.c = trial_data$comparator_cs,
      n.c = trial_data$comparator_n,
      studlab = trial_data$label,
      sm = "RR",
      method = "Inverse",
      method.tau = "DL",
      method.random.ci = "HK",
      incr = continuity_correction,
      method.incr = "only0",
      random = TRUE,
      common = FALSE,
      subgroup = trial_data$family_label
    )
  )
  log_line(
    sprintf(
      paste0(
        "  incr = %.1f : overall RR = %.3f (%.3f to %.3f); ",
        "family A RR = %.3f (%.3f to %.3f)"
      ),
      continuity_correction,
      exp(continuity_model$TE.random),
      exp(continuity_model$lower.random),
      exp(continuity_model$upper.random),
      exp(continuity_model$TE.random.w[1]),
      exp(continuity_model$lower.random.w[1]),
      exp(continuity_model$upper.random.w[1])
    )
  )
}
