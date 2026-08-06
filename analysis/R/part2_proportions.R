###############################################################################
# PART 2. SINGLE-ARM POOLED PROPORTIONS BY SCREENING STRATEGY
###############################################################################
log_section("PART 2. POOLED PROPORTIONS BY SCREENING STRATEGY (GLMM)")
log_line(
  "Arms are classified by the triage test that triggers work-up. Denominators"
)
log_line("are the number of men actually screened where reported; control arms of")
log_line("invitation trials were not screened and retain the invited/randomized N.")
log_line("")
log_line("IMPORTANT CAVEAT: in ERSPC, PLCO and CAP the numerators are cancers")
log_line("accrued across the whole randomized arm over 9-10 years of follow-up,")
log_line("whereas the per-screened denominator counts only men actually tested.")
log_line("These per-screened rates are therefore upward-biased for those three")
log_line("trials. Part 3 repeats the analysis on ITT denominators.")
log_line("")
log_line("Arms per strategy group:")

for (strategy in levels(droplevels(arm_data$strategy))) {
  uses_strategy <- arm_data$strategy == strategy
  arm_labels <- paste0(
    arm_data$trial[uses_strategy],
    "-",
    arm_data$arm[uses_strategy]
  )

  log_line(
    sprintf(
      "  %-28s %s",
      strategy,
      paste(arm_labels, collapse = ", ")
    )
  )
}

gg1_per_screened_proportion <- fit_proportion(
  source_data = arm_data,
  event_column = "gg1",
  total_column = "n_screened",
  title = "ISUP GG1 detection rate per man screened, by screening strategy"
)
report_proportion(
  gg1_per_screened_proportion,
  "2.1 PRIMARY - GG1 detection rate per man screened, by strategy"
)
report_leave_one_out_proportion(
  gg1_per_screened_proportion,
  "GG1 rate per screened man"
)
collect_proportion_results(
  gg1_per_screened_proportion,
  "GG1 rate per man screened (proportion)"
)
save_proportion_forest_plot(
  gg1_per_screened_proportion,
  "fig6_prop_gg1_per_screened",
  "ISUP GG1 detected per 100 men screened (%)"
)

cspca_per_screened_proportion <- fit_proportion(
  source_data = arm_data,
  event_column = "cs",
  total_column = "n_screened",
  title = "ISUP >=2 detection rate per man screened, by screening strategy"
)
report_proportion(
  cspca_per_screened_proportion,
  "2.2 SECONDARY - csPCa detection rate per man screened, by strategy"
)
report_leave_one_out_proportion(
  cspca_per_screened_proportion,
  "csPCa rate per screened man"
)
collect_proportion_results(
  cspca_per_screened_proportion,
  "csPCa rate per man screened (proportion)"
)
save_proportion_forest_plot(
  cspca_per_screened_proportion,
  "fig7_prop_cspca_per_screened",
  "ISUP >=2 detected per 100 men screened (%)"
)

gg1_share_proportion <- fit_proportion(
  source_data = arm_data,
  event_column = "gg1",
  total_column = "pca",
  title = paste0(
    "ISUP GG1 as a proportion of all detected ",
    "prostate cancers, by strategy"
  )
)
report_proportion(
  gg1_share_proportion,
  "2.3 PRIMARY (grade mix) - GG1 share of all detected PCa, by strategy"
)
report_leave_one_out_proportion(
  gg1_share_proportion,
  "GG1 share of detected PCa"
)
collect_proportion_results(
  gg1_share_proportion,
  "GG1 share of detected PCa (proportion)"
)
save_proportion_forest_plot(
  gg1_share_proportion,
  "fig8_prop_gg1_share_of_pca",
  "ISUP GG1 as % of all detected prostate cancers"
)

biopsies_per_screened_proportion <- fit_proportion(
  source_data = arm_data,
  event_column = "bx",
  total_column = "n_screened",
  title = "Prostate biopsies per man screened, by screening strategy"
)
report_proportion(
  biopsies_per_screened_proportion,
  "2.4 SECONDARY - Biopsies per man screened, by strategy"
)
log_line(
  paste0(
    "  Arms with biopsy counts NR are excluded ",
    "(ERSPC-C, PLCO both arms, CAP-C, ProScreen-C)."
  )
)
report_leave_one_out_proportion(
  biopsies_per_screened_proportion,
  "biopsies per screened man"
)
collect_proportion_results(
  biopsies_per_screened_proportion,
  "Biopsies per man screened (proportion)"
)
save_proportion_forest_plot(
  biopsies_per_screened_proportion,
  "fig9_prop_biopsies_per_screened",
  "Biopsies per 100 men screened (%)"
)
