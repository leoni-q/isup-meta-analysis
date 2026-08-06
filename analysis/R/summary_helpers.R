# Shared result formatting and summary-table collection.

# Result formatting ---------------------------------------------------------

format_heterogeneity <- function(model) {
  # GLMM fits return Q / df.Q / pval.Q of length 2: [1] Cochran's Q (Wald-type),
  # [2] the likelihood-ratio test of heterogeneity. Protocol specifies Cochran's Q.
  cochrans_q <- model$Q[1]
  q_degrees_of_freedom <- model$df.Q[1]
  q_p_value <- model$pval.Q[1]

  q_degrees_of_freedom_text <- if (is.na(q_degrees_of_freedom)) {
    "NA"
  } else {
    as.character(round(q_degrees_of_freedom))
  }
  q_p_value_text <- if (is.na(q_p_value)) {
    "NA"
  } else {
    format.pval(q_p_value, digits = 4, eps = 1e-4)
  }
  q_summary <- sprintf(
    "Cochran's Q = %.2f, df = %s, p = %s",
    cochrans_q,
    q_degrees_of_freedom_text,
    q_p_value_text
  )

  i_squared_summary <- if (is.na(model$I2[1])) {
    "I2 = NA"
  } else {
    sprintf("I2 = %.1f%%", 100 * model$I2[1])
  }
  tau_squared_summary <- if (is.na(model$tau2[1])) {
    "tau2 = NA"
  } else {
    sprintf("tau2 = %.4f", model$tau2[1])
  }
  significance_note <- if (!is.na(q_p_value) && q_p_value < 0.05) {
    "  ** significant heterogeneity **"
  } else {
    ""
  }
  likelihood_ratio_note <- if (
    length(model$pval.Q) > 1 && !is.na(model$pval.Q[2])
  ) {
    sprintf(
      "; LRT of heterogeneity p = %s",
      format.pval(model$pval.Q[2], digits = 4, eps = 1e-4)
    )
  } else {
    ""
  }

  paste0(
    "  Heterogeneity: ",
    q_summary,
    "; ",
    i_squared_summary,
    "; ",
    tau_squared_summary,
    likelihood_ratio_note,
    significance_note
  )
}

first_or_na <- function(values) {
  if (length(values) == 0) {
    return(NA)
  }

  values[1]
}

# Summary-table collection ---------------------------------------------------

# Summary rows are accumulated as each model runs. `<<-` assigns into the
# enclosing script scope (the R equivalent of mutating a captured variable).
summary_rows <- list()

append_summary_row <- function(
  analysis,
  scale,
  group,
  study_count,
  estimate,
  confidence_interval_lower,
  confidence_interval_upper,
  cochrans_q,
  q_degrees_of_freedom,
  q_p_value,
  i_squared_percent
) {
  summary_rows[[length(summary_rows) + 1]] <<- data.frame(
    analysis = analysis,
    scale = scale,
    group = group,
    k = study_count,
    estimate = estimate,
    ci_lower = confidence_interval_lower,
    ci_upper = confidence_interval_upper,
    Q = cochrans_q,
    df_Q = q_degrees_of_freedom,
    p_Q = q_p_value,
    I2_pct = i_squared_percent,
    stringsAsFactors = FALSE
  )
}

collect_risk_ratio_results <- function(model, analysis_name) {
  append_summary_row(
    analysis_name,
    "RR",
    "Overall (not interpreted across families)",
    model$k,
    exp(model$TE.random),
    exp(model$lower.random),
    exp(model$upper.random),
    model$Q,
    model$df.Q,
    model$pval.Q,
    100 * model$I2
  )

  if (!is.null(model$subgroup)) {
    for (subgroup_index in seq_along(model$subgroup.levels)) {
      append_summary_row(
        analysis_name,
        "RR",
        model$subgroup.levels[subgroup_index],
        model$k.w[subgroup_index],
        exp(model$TE.random.w[subgroup_index]),
        exp(model$lower.random.w[subgroup_index]),
        exp(model$upper.random.w[subgroup_index]),
        model$Q.w[subgroup_index],
        model$df.Q.w[subgroup_index],
        model$pval.Q.w[subgroup_index],
        100 * model$I2.w[subgroup_index]
      )
    }
  }
}

collect_proportion_results <- function(model, analysis_name) {
  append_summary_row(
    analysis_name,
    "Proportion (%)",
    "All arms pooled",
    model$k,
    100 * plogis(model$TE.random),
    100 * plogis(model$lower.random),
    100 * plogis(model$upper.random),
    first_or_na(model$Q),
    first_or_na(model$df.Q),
    first_or_na(model$pval.Q),
    100 * first_or_na(model$I2)
  )

  if (!is.null(model$subgroup)) {
    for (subgroup_index in seq_along(model$subgroup.levels)) {
      if (model$k.w[subgroup_index] == 1) {
        arm_index <- which(
          as.character(model$subgroup) == model$subgroup.levels[subgroup_index]
        )
        exact_confidence_interval <- binom.test(
          model$event[arm_index],
          model$n[arm_index]
        )$conf.int

        append_summary_row(
          analysis_name,
          "Proportion (%)",
          model$subgroup.levels[subgroup_index],
          1,
          100 * model$event[arm_index] / model$n[arm_index],
          100 * exact_confidence_interval[1],
          100 * exact_confidence_interval[2],
          NA,
          NA,
          NA,
          NA
        )
      } else {
        append_summary_row(
          analysis_name,
          "Proportion (%)",
          model$subgroup.levels[subgroup_index],
          model$k.w[subgroup_index],
          100 * plogis(model$TE.random.w[subgroup_index]),
          100 * plogis(model$lower.random.w[subgroup_index]),
          100 * plogis(model$upper.random.w[subgroup_index]),
          first_or_na(model$Q.w[subgroup_index]),
          first_or_na(model$df.Q.w[subgroup_index]),
          first_or_na(model$pval.Q.w[subgroup_index]),
          100 * first_or_na(model$I2.w[subgroup_index])
        )
      }
    }
  }
}
