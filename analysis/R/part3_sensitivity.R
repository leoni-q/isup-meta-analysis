###############################################################################
# PART 3. SENSITIVITY ANALYSES
###############################################################################
log_section("PART 3. SENSITIVITY ANALYSES")
sensitivity_models <- list()

run_rr_sensitivity_set <- function(
  source_data,
  sensitivity_name,
  group_name,
  outcomes = RR_OUTCOMES
) {
  fitted_models <- list()

  for (outcome_name in names(outcomes)) {
    outcome <- outcomes[[outcome_name]]
    complete_data <- complete_rr_data(
      source_data,
      outcome$intervention_event_column,
      outcome$intervention_total_column,
      outcome$comparator_event_column,
      outcome$comparator_total_column
    )
    if (nrow(complete_data) == 0) {
      log_line(
        sprintf(
          "  Skipping %s: no trial reports both arms.",
          outcome$analysis_name
        )
      )
      next
    }

    model <- fit_risk_ratio(
      source_data = source_data,
      intervention_event_column = outcome$intervention_event_column,
      intervention_total_column = outcome$intervention_total_column,
      comparator_event_column = outcome$comparator_event_column,
      comparator_total_column = outcome$comparator_total_column,
      title = paste(outcome$analysis_name, sensitivity_name, sep = " - "),
      include_subgroups = FALSE
    )
    fitted_models[[outcome_name]] <- model
    report_risk_ratio(
      model,
      paste(sensitivity_name, "-", toupper(outcome$analysis_name))
    )
    collect_risk_ratio_results(
      model,
      paste0("SENS ", sensitivity_name, ": ", outcome$analysis_name),
      group_name
    )
  }

  fitted_models
}

# 3.1 CAP cluster-randomized design: unadjusted CAP remains in the primary
# invitation analysis, so a full excluding-CAP sensitivity is reported.
log_line("")
log_line("--- 3.1 Invitation-family RRs excluding CAP ---")
invitation_without_cap <- invitation_trial_data[
  invitation_trial_data$trial != "CAP",
]
sensitivity_models$without_cap <- run_rr_sensitivity_set(
  source_data = invitation_without_cap,
  sensitivity_name = "excluding CAP",
  group_name = paste(FAMILY_LABELS[["A"]], "(CAP excluded)")
)

# 3.2 ERSPC Gleason counts were reconstructed from percentages.
log_line("")
log_line("--- 3.2 Invitation-family RRs excluding ERSPC ---")
invitation_without_erspc <- invitation_trial_data[
  invitation_trial_data$trial != "ERSPC",
]
sensitivity_models$without_erspc <- run_rr_sensitivity_set(
  source_data = invitation_without_erspc,
  sensitivity_name = "excluding ERSPC",
  group_name = paste(FAMILY_LABELS[["A"]], "(ERSPC excluded)")
)

# 3.3 Main diagnostic-strategy RRs use randomized denominators. This sensitivity
# uses the assessed/screened denominators for MVP and PROSA and is therefore a
# genuine alternative rather than the previous no-op reassignment.
log_line("")
log_line(
  "--- 3.3 Diagnostic-family detection RRs using assessed denominators ---"
)
diagnostic_assessed_data <- diagnostic_trial_data
uses_assessed_denominator <- diagnostic_assessed_data$trial %in% c("MVP", "PROSA")
diagnostic_assessed_data$intervention_n[uses_assessed_denominator] <- (
  diagnostic_assessed_data$intervention_n_screened[uses_assessed_denominator]
)
diagnostic_assessed_data$comparator_n[uses_assessed_denominator] <- (
  diagnostic_assessed_data$comparator_n_screened[uses_assessed_denominator]
)

detection_outcomes <- RR_OUTCOMES[c("gg1_detection", "gg2plus_detection")]
sensitivity_models$assessed_denominators <- run_rr_sensitivity_set(
  source_data = diagnostic_assessed_data,
  sensitivity_name = "MVP/PROSA assessed denominators",
  group_name = paste(FAMILY_LABELS[["B"]], "(assessed-denominator sensitivity)"),
  outcomes = detection_outcomes
)

# 3.4 PROBASE is the only retained contrast with a zero GG>=2 control cell.
# Vary the all-cell correction for this standalone timing comparison.
log_line("")
log_line(
  "--- 3.4 PROBASE GG>=2 continuity-correction sensitivity (all-cell rule) ---"
)

for (continuity_correction in c(0.1, 0.5, 1)) {
  continuity_model <- suppressWarnings(
    metabin(
      event.e = timing_trial_data$intervention_cs,
      n.e = timing_trial_data$intervention_n,
      event.c = timing_trial_data$comparator_cs,
      n.c = timing_trial_data$comparator_n,
      studlab = timing_trial_data$label,
      sm = "RR",
      method = "Inverse",
      method.tau = "DL",
      method.random.ci = "HK",
      incr = continuity_correction,
      method.incr = "only0",
      RR.Cochrane = TRUE,
      random = TRUE,
      common = FALSE
    )
  )
  log_line(
    sprintf(
      "  incr = %.1f : PROBASE RR = %.3f (%.3f to %.3f)",
      continuity_correction,
      exp(continuity_model$TE.random),
      exp(continuity_model$lower.random),
      exp(continuity_model$upper.random)
    )
  )
}
