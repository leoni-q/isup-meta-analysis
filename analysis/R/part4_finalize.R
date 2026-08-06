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
    "comparative analysis here contains k = %d trials, and every strategy",
    nrow(trial_data)
  )
)
log_line(
  "subgroup contains fewer still. Publication bias was therefore NOT assessed,"
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
write.csv(
  summary_table,
  file.path(OUTPUT_DIR, "table2_pooled_estimates.csv"),
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
