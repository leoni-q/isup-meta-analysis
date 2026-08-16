# Structural validation for generated meta-analysis outputs.

command_arguments <- commandArgs(trailingOnly = FALSE)
script_argument <- command_arguments[grepl("^--file=", command_arguments)]
script_path <- sub("^--file=", "", script_argument[1])
ANALYSIS_DIR <- dirname(normalizePath(script_path))
OUTPUT_DIR <- file.path(ANALYSIS_DIR, "output")

arm_data <- read.csv(
  file.path(ANALYSIS_DIR, "arm_level_data.csv"),
  stringsAsFactors = FALSE
)
summary_table <- read.csv(
  file.path(OUTPUT_DIR, "table2_pooled_estimates.csv"),
  stringsAsFactors = FALSE
)
leave_one_out_table <- read.csv(
  file.path(OUTPUT_DIR, "table3_leave_one_out.csv"),
  stringsAsFactors = FALSE
)
study_strategy_units <- read.csv(
  file.path(OUTPUT_DIR, "table1b_study_strategy_units.csv"),
  stringsAsFactors = FALSE
)

# Cohort membership and PROBASE endpoint classification.
stopifnot(
  all(arm_data$family[arm_data$trial == "PROBASE"] == "T"),
  all(arm_data$family[arm_data$trial == "CAP"] == "A"),
  is.na(arm_data$n_screened[
    arm_data$trial == "PROBASE" & arm_data$arm == "C"
  ])
)

# CAP uses the mature 15-year Martin 2024 extraction from the supplied
# primary-outcome CSV, rather than the superseded 10-year Martin 2018 counts.
cap_data <- arm_data[arm_data$trial == "CAP", ]
cap_data <- cap_data[match(c("I", "C"), cap_data$arm), ]
stopifnot(
  identical(cap_data$year, c(2024L, 2024L)),
  identical(cap_data$n_itt, c(189326L, 219395L)),
  identical(cap_data$gg1, c(4111L, 3482L)),
  identical(cap_data$cs, c(6649L, 8019L)),
  identical(cap_data$pca, c(12013L, 12958L))
)

# One independent contribution per trial and strategy.
stopifnot(
  !any(duplicated(study_strategy_units[c("trial", "strategy")])),
  nrow(study_strategy_units[study_strategy_units$trial == "PROBASE", ]) == 1,
  study_strategy_units$source_arms[
    study_strategy_units$trial == "PROBASE"
  ] == "I",
  study_strategy_units$arm_count[
    study_strategy_units$trial == "GOTEBORG-2"
  ] == 2
)

# Detection-rate numerators are nested in whole randomized/invited arms.
stopifnot(
  all(study_strategy_units$gg1 <= study_strategy_units$n_itt),
  all(study_strategy_units$cs <= study_strategy_units$n_itt),
  all(study_strategy_units$gg1 + study_strategy_units$cs <=
    study_strategy_units$pca),
  all(study_strategy_units$gg1[
    !is.na(study_strategy_units$gg1_bx_n)
  ] <= study_strategy_units$gg1_bx_n[
    !is.na(study_strategy_units$gg1_bx_n)
  ]),
  all(study_strategy_units$cs[
    !is.na(study_strategy_units$cs_bx_n)
  ] <= study_strategy_units$cs_bx_n[
    !is.na(study_strategy_units$cs_bx_n)
  ])
)

# Biopsy-yield denominators follow the supplied outcome cells rather than being
# inferred from bx. This excludes CAP's incompatible first-round biopsy count
# and preserves the supplied STHLM3-MRI GG>=2 denominator of 331.
stopifnot(
  is.na(arm_data$gg1_bx_n[
    arm_data$trial == "CAP" & arm_data$arm == "I"
  ]),
  is.na(arm_data$cs_bx_n[
    arm_data$trial == "CAP" & arm_data$arm == "I"
  ]),
  arm_data$gg1_bx_n[
    arm_data$trial == "STHLM3-MRI" & arm_data$arm == "I"
  ] == 311,
  arm_data$cs_bx_n[
    arm_data$trial == "STHLM3-MRI" & arm_data$arm == "I"
  ] == 331
)

# Main RR models contain the six agreed outcomes and no cross-family pooled
# estimate. Biopsy-yield RRs are absent from the invitation family because no
# invitation trial reports compatible biopsy denominators for both arms.
main_rr <- summary_table[
  summary_table$scale == "RR" &
    !grepl("^SENS ", summary_table$analysis),
]
stopifnot(
  identical(
    sort(unique(main_rr$analysis)),
    sort(
      c(
        "GG1 detection per randomized/invited participant",
        "GG1 among all detected prostate cancers",
        "GG>=2 detection per randomized/invited participant",
        "GG>=2 among all detected prostate cancers",
        "GG1 among participants biopsied",
        "GG>=2 among participants biopsied"
      )
    )
  ),
  !any(grepl("Overall", main_rr$group, fixed = TRUE)),
  all(main_rr$k[
    grepl("^A\\.", main_rr$group)
  ] == 4),
  all(main_rr$k[
    grepl("^B\\.", main_rr$group)
  ] == 4),
  all(main_rr$k[
    grepl("^T\\.", main_rr$group)
  ] == 1),
  !any(
    grepl("participants biopsied", main_rr$analysis, fixed = TRUE) &
      grepl("^A\\.", main_rr$group)
  )
)

# Every pooled estimate with at least two studies omits every contributing
# study exactly once. Single-study estimates correctly have no LOO row.
main_summary <- summary_table[!grepl("^SENS ", summary_table$analysis), ]
for (row_index in seq_len(nrow(main_summary))) {
  result_row <- main_summary[row_index, ]
  matching_loo <- leave_one_out_table[
    leave_one_out_table$analysis == result_row$analysis &
      leave_one_out_table$group == result_row$group,
  ]
  expected_rows <- if (result_row$k >= 2) result_row$k else 0
  stopifnot(
    nrow(matching_loo) == expected_rows,
    !any(duplicated(matching_loo$omitted_study))
  )
}

# Cochran Q degrees of freedom are group-specific, not reused overall values.
pooled_rows <- summary_table$k >= 2
stopifnot(
  all(summary_table$df_Q[pooled_rows] == summary_table$k[pooled_rows] - 1)
)

message(
  sprintf(
    "Validation passed: %d pooled rows, %d leave-one-out refits, %d study-strategy units.",
    nrow(summary_table),
    nrow(leave_one_out_table),
    nrow(study_strategy_units)
  )
)
