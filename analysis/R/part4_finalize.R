###############################################################################
# PART 4. PUBLICATION BIAS
###############################################################################
log_section("PART 4. PUBLICATION BIAS")
log_line(
  "The protocol restricts publication-bias assessment (Peters' regression and"
)
log_line("funnel plots) to analyses containing 10 or more studies. The largest")
log_line(
  sprintf(
    "reported pooled analysis here contains k = %d studies, and every model",
    max(
      c(
        vapply(risk_ratio_models, function(model) model$k, numeric(1)),
        vapply(proportion_models, function(model) model$k, numeric(1))
      )
    )
  )
)
log_line(
  "contains fewer than 10 studies. Publication bias was therefore NOT assessed,"
)
log_line("and small-study effects cannot be excluded.")

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
