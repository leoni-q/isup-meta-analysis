###############################################################################
# PART 1. COMPARATIVE ANALYSES - RISK RATIOS (ITT denominators)
###############################################################################
log_section(
  "PART 1. COMPARATIVE META-ANALYSES (RISK RATIOS, ITT DENOMINATORS)"
)
log_line("Denominators are randomized / invited N throughout, preserving the")
log_line(
  "randomized comparison. Subgrouped by contrast family; the overall pooled"
)
log_line(
  "RR across both families is computed but NOT interpreted (protocol s.8/s.9)."
)

gg1_detection_risk_ratio <- fit_risk_ratio(
  source_data = trial_data,
  intervention_event_column = "intervention_gg1",
  intervention_total_column = "intervention_n",
  comparator_event_column = "comparator_gg1",
  comparator_total_column = "comparator_n",
  title = paste0(
    "Primary outcome: ISUP grade group 1 detection ",
    "per man randomized/invited"
  )
)
report_risk_ratio(
  gg1_detection_risk_ratio,
  "1.1 PRIMARY OUTCOME - ISUP GG1 detection per man randomized/invited (RR)"
)
report_leave_one_out_risk_ratio(
  gg1_detection_risk_ratio,
  "GG1 detection RR"
)
collect_risk_ratio_results(
  gg1_detection_risk_ratio,
  "GG1 detection per man (RR)"
)
save_risk_ratio_forest_plot(
  gg1_detection_risk_ratio,
  "fig1_rr_gg1_detection",
  "Favours intervention",
  "Favours comparator"
)

cspca_detection_risk_ratio <- fit_risk_ratio(
  source_data = trial_data,
  intervention_event_column = "intervention_cs",
  intervention_total_column = "intervention_n",
  comparator_event_column = "comparator_cs",
  comparator_total_column = "comparator_n",
  title = paste0(
    "Secondary outcome: ISUP grade group >=2 detection ",
    "per man randomized/invited"
  )
)
report_risk_ratio(
  cspca_detection_risk_ratio,
  "1.2 SECONDARY - ISUP >=2 (csPCa) detection per man randomized/invited (RR)"
)
report_leave_one_out_risk_ratio(
  cspca_detection_risk_ratio,
  "csPCa detection RR"
)
collect_risk_ratio_results(
  cspca_detection_risk_ratio,
  "csPCa detection per man (RR)"
)
save_risk_ratio_forest_plot(
  cspca_detection_risk_ratio,
  "fig2_rr_cspca_detection",
  "Fewer csPCa",
  "More csPCa",
  x_ticks = c(0.5, 1, 2, 5, 10, 20, 50, 100)
)

any_pca_detection_risk_ratio <- fit_risk_ratio(
  source_data = trial_data,
  intervention_event_column = "intervention_pca",
  intervention_total_column = "intervention_n",
  comparator_event_column = "comparator_pca",
  comparator_total_column = "comparator_n",
  title = "Any prostate cancer detection per man randomized/invited"
)
report_risk_ratio(
  any_pca_detection_risk_ratio,
  "1.3 CONTEXT - Any prostate cancer detection per man randomized/invited (RR)"
)
report_leave_one_out_risk_ratio(
  any_pca_detection_risk_ratio,
  "any PCa detection RR"
)
collect_risk_ratio_results(
  any_pca_detection_risk_ratio,
  "Any PCa detection per man (RR)"
)
save_risk_ratio_forest_plot(
  any_pca_detection_risk_ratio,
  "fig3_rr_any_pca",
  "Fewer PCa",
  "More PCa",
  x_ticks = c(0.5, 1, 2, 5, 10, 25, 50)
)

gg1_share_risk_ratio <- fit_risk_ratio(
  source_data = trial_data,
  intervention_event_column = "intervention_gg1",
  intervention_total_column = "intervention_pca",
  comparator_event_column = "comparator_gg1",
  comparator_total_column = "comparator_pca",
  title = paste0(
    "Grade mix: ISUP GG1 as a proportion of ",
    "all detected prostate cancers"
  )
)
report_risk_ratio(
  gg1_share_risk_ratio,
  "1.4 PRIMARY (grade mix) - GG1 as a share of all detected PCa (RR)"
)
report_leave_one_out_risk_ratio(
  gg1_share_risk_ratio,
  "GG1 share of detected PCa RR"
)
collect_risk_ratio_results(
  gg1_share_risk_ratio,
  "GG1 share of detected PCa (RR)"
)
save_risk_ratio_forest_plot(
  gg1_share_risk_ratio,
  "fig4_rr_gg1_share_of_pca",
  "Lower GG1 share",
  "Higher GG1 share",
  x_ticks = c(0.2, 0.35, 0.5, 0.7, 1, 1.5, 2)
)

biopsy_risk_ratio <- fit_risk_ratio(
  source_data = trial_data,
  intervention_event_column = "intervention_biopsies",
  intervention_total_column = "intervention_n",
  comparator_event_column = "comparator_biopsies",
  comparator_total_column = "comparator_n",
  title = "Secondary outcome: prostate biopsies performed per man randomized"
)
report_risk_ratio(
  biopsy_risk_ratio,
  "1.5 SECONDARY - Biopsies performed per man randomized (RR)"
)
log_line("  NOTE: restricted to trials reporting biopsies in BOTH arms.")
log_line(
  paste0(
    "  ERSPC, CAP and ProScreen report intervention-arm biopsies only; ",
    "PLCO reports neither."
  )
)
report_leave_one_out_risk_ratio(biopsy_risk_ratio, "biopsy RR")
collect_risk_ratio_results(biopsy_risk_ratio, "Biopsies per man (RR)")
save_risk_ratio_forest_plot(
  biopsy_risk_ratio,
  "fig5_rr_biopsies",
  "Fewer biopsies",
  "More biopsies",
  x_ticks = c(0.25, 0.5, 1, 2, 4)
)
