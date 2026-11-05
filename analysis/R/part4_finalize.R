###############################################################################
# PART 4. PUBLICATION BIAS
###############################################################################
log_section("PART 4. PUBLICATION BIAS")
log_line(
  "The protocol restricts publication-bias assessment (Peters' regression and"
)
log_line("funnel plots) to analyses containing 10 or more studies.")

PUBLICATION_BIAS_MINIMUM_STUDIES <- 10
all_models <- c(risk_ratio_models, proportion_models)
model_study_counts <- vapply(all_models, function(model) model$k, numeric(1))
eligible_models <- all_models[
  model_study_counts >= PUBLICATION_BIAS_MINIMUM_STUDIES
]

log_line(
  sprintf(
    "The largest reported pooled analysis contains k = %d studies.",
    max(model_study_counts)
  )
)

# The threshold is evaluated rather than asserted, so adding trials in a search
# update triggers the mandated test instead of silently keeping this statement.
if (length(eligible_models) == 0) {
  log_line(
    sprintf(
      "No model reaches k >= %d, so publication bias was NOT assessed and",
      PUBLICATION_BIAS_MINIMUM_STUDIES
    )
  )
  log_line("small-study effects cannot be excluded.")
} else {
  for (model_name in names(eligible_models)) {
    model <- eligible_models[[model_name]]
    peters_test <- suppressWarnings(
      try(metabias(model, method.bias = "Peters"), silent = TRUE)
    )
    if (inherits(peters_test, "try-error")) {
      log_line(sprintf("  %s: Peters' test could not be computed.", model_name))
      next
    }
    log_line(
      sprintf(
        "  %s: Peters' test t = %.3f, p = %s (k = %d)",
        model_name,
        peters_test$statistic,
        format.pval(peters_test$pval, digits = 4, eps = 1e-4),
        model$k
      )
    )
    funnel_path <- file.path(
      OUTPUT_DIR,
      paste0("funnel_", model_name, ".png")
    )
    png(funnel_path, width = 7, height = 6, units = "in", res = 200)
    funnel(model)
    dev.off()
  }
}

###############################################################################
# Protocol conformance
###############################################################################
log_section("PROTOCOL DEVIATIONS AND UNRESOLVED ITEMS")
log_line(
  "Deviations from protocol.md are recorded here so they can be reported"
)
log_line("alongside the results rather than discovered on comparison.")
log_line("")
log_line("1. Stratification of comparative analyses (protocol s.8).")
log_line("   The protocol specifies pooling risk ratios 'stratified by the")
log_line("   screening strategy'. A risk ratio is a contrast BETWEEN two")
log_line("   strategies, so a single strategy label cannot describe it. Trials")
log_line("   are therefore grouped by contrast family: screening invitation vs")
log_line("   no screening (A), diagnostic/biopsy strategy with both arms")
log_line("   screened (B), and timing of screening initiation (T). Families are")
log_line("   never pooled with one another. Single-arm proportions ARE")
log_line("   stratified by screening strategy exactly as the protocol states.")
log_line("")
log_line("2. Continuity correction for proportions (protocol s.8).")
log_line("   The protocol states a 0.5 correction for zero cells in both")
log_line("   analyses. It is applied to the risk-ratio models. The proportion")
log_line("   models are logit-normal GLMMs fitted by likelihood, which use the")
log_line("   exact binomial likelihood and require no continuity correction;")
log_line("   adding one would bias them. No correction is therefore applied on")
log_line("   the proportion side.")
log_line("")
log_line("3. Curative-treatment outcome (protocol s.2 and s.9).")
log_line("   'Percentage of patients receiving curative treatment over the")
log_line("   course of follow-up' is a registered secondary outcome. The")
log_line("   included reports do not provide extractable arm-level treatment")
log_line("   counts, so no data column exists and the outcome is NOT analysed.")
log_line("   This is an outcome-level deviation and must be reported as such.")
log_line("")
log_line("4. Adverse events and procedural harms (protocol s.9).")
log_line("   Contextual data only, not formal endpoints; not extracted into the")
log_line("   analysis dataset and therefore not summarized here.")
log_line("")
log_line("5. Exclusion of an extracted trial.")
log_line("   ProScreen was extracted (arm_level_data.csv) but is excluded from")
log_line("   every analysis, including the descriptive tables.")

###############################################################################
# Export summary table
###############################################################################
summary_table <- do.call(rbind, summary_rows)
summary_table$estimate <- round(summary_table$estimate, 4)
summary_table$ci_lower <- round(summary_table$ci_lower, 4)
summary_table$ci_upper <- round(summary_table$ci_upper, 4)
summary_table$Q <- round(summary_table$Q, 3)
summary_table$p_Q <- signif(summary_table$p_Q, 4)
summary_table$I2_pct <- round(summary_table$I2_pct, 1)
summary_table$pi_lower <- round(summary_table$pi_lower, 4)
summary_table$pi_upper <- round(summary_table$pi_upper, 4)
write.csv(
  summary_table,
  file.path(OUTPUT_DIR, "table2_pooled_estimates.csv"),
  row.names = FALSE
)

leave_one_out_table <- do.call(rbind, leave_one_out_rows)
leave_one_out_table$estimate <- round(leave_one_out_table$estimate, 4)
leave_one_out_table$ci_lower <- round(leave_one_out_table$ci_lower, 4)
leave_one_out_table$ci_upper <- round(leave_one_out_table$ci_upper, 4)
leave_one_out_table$I2_pct <- round(leave_one_out_table$I2_pct, 1)
write.csv(
  leave_one_out_table,
  file.path(OUTPUT_DIR, "table3_leave_one_out.csv"),
  row.names = FALSE
)

write.csv(
  proportion_data,
  file.path(OUTPUT_DIR, "table1b_study_strategy_units.csv"),
  row.names = FALSE
)

log_section("SESSION INFO")
log_line(
  sprintf(
    "R %s.%s | meta %s | metafor %s",
    R.version$major,
    R.version$minor,
    packageVersion("meta"),
    packageVersion("metafor")
  )
)
log_line(sprintf("Outputs written to %s", OUTPUT_DIR))
