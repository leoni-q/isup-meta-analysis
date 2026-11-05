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

format_forest_heterogeneity_label <- function(model) {
  # Heterogeneity line drawn under each forest plot. meta::forest() can only
  # label the Q statistic as chi-squared, so the plots use this custom text:
  # Cochran's Q test (statistic, df, p) plus tau^2, with I^2 deliberately
  # omitted. GLMM fits return Q of length 2; [1] is Cochran's Q.
  cochrans_q <- model$Q[1]
  q_degrees_of_freedom <- model$df.Q[1]
  q_p_value <- model$pval.Q[1]

  if (is.na(cochrans_q) || is.na(q_degrees_of_freedom)) {
    return("Heterogeneity: Cochran's Q not estimable (k = 1)")
  }

  p_value_text <- if (is.na(q_p_value)) {
    "p = NA"
  } else if (q_p_value < 0.0001) {
    "p < 0.0001"
  } else {
    sprintf("p = %.4f", q_p_value)
  }
  tau_squared_text <- if (is.na(model$tau2[1])) {
    ""
  } else {
    sprintf("; Tau\u00b2 = %.4f", model$tau2[1])
  }

  sprintf(
    "Heterogeneity: Cochran's Q = %.2f, df = %d, %s%s",
    cochrans_q,
    as.integer(round(q_degrees_of_freedom)),
    p_value_text,
    tau_squared_text
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
leave_one_out_rows <- list()

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
  i_squared_percent,
  prediction_interval_lower = NA_real_,
  prediction_interval_upper = NA_real_
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
    pi_lower = prediction_interval_lower,
    pi_upper = prediction_interval_upper,
    stringsAsFactors = FALSE
  )
}

append_leave_one_out_row <- function(
  analysis,
  group,
  omitted_study,
  remaining_studies,
  scale,
  estimate,
  confidence_interval_lower,
  confidence_interval_upper,
  i_squared_percent = NA_real_
) {
  leave_one_out_rows[[length(leave_one_out_rows) + 1]] <<- data.frame(
    analysis = analysis,
    group = group,
    omitted_study = omitted_study,
    k_remaining = remaining_studies,
    scale = scale,
    estimate = estimate,
    ci_lower = confidence_interval_lower,
    ci_upper = confidence_interval_upper,
    I2_pct = i_squared_percent,
    stringsAsFactors = FALSE
  )
}

collect_risk_ratio_results <- function(model, analysis_name, group_name) {
  prediction_interval_lower <- if (
    model$k >= 3 && length(model$lower.predict) > 0
  ) {
    exp(first_or_na(model$lower.predict))
  } else {
    NA_real_
  }
  prediction_interval_upper <- if (
    model$k >= 3 && length(model$upper.predict) > 0
  ) {
    exp(first_or_na(model$upper.predict))
  } else {
    NA_real_
  }

  # metabin returns Q = 0, df = 0 for a single contrast, which would read as a
  # heterogeneity test that was actually never run. The proportion collector
  # already reports NA in that situation; match it so the exported table
  # encodes single-study rows the same way on both scales.
  has_heterogeneity_test <- model$k > 1
  cochrans_q <- if (has_heterogeneity_test) model$Q else NA_real_
  q_degrees_of_freedom <- if (has_heterogeneity_test) model$df.Q else NA_real_
  q_p_value <- if (has_heterogeneity_test) model$pval.Q else NA_real_
  i_squared_percent <- if (has_heterogeneity_test) 100 * model$I2 else NA_real_

  append_summary_row(
    analysis_name,
    "RR",
    group_name,
    model$k,
    exp(model$TE.random),
    exp(model$lower.random),
    exp(model$upper.random),
    cochrans_q,
    q_degrees_of_freedom,
    q_p_value,
    i_squared_percent,
    prediction_interval_lower,
    prediction_interval_upper
  )
}

extract_proportion_result <- function(model) {
  if (model$k == 1) {
    exact_confidence_interval <- binom.test(model$event[1], model$n[1])$conf.int
    return(
      list(
        estimate = 100 * model$event[1] / model$n[1],
        lower = 100 * exact_confidence_interval[1],
        upper = 100 * exact_confidence_interval[2],
        Q = NA_real_,
        df_Q = NA_real_,
        p_Q = NA_real_,
        I2 = NA_real_,
        pi_lower = NA_real_,
        pi_upper = NA_real_
      )
    )
  }

  list(
    estimate = 100 * plogis(model$TE.random),
    lower = 100 * plogis(model$lower.random),
    upper = 100 * plogis(model$upper.random),
    Q = first_or_na(model$Q),
    df_Q = first_or_na(model$df.Q),
    p_Q = first_or_na(model$pval.Q),
    I2 = 100 * first_or_na(model$I2),
    pi_lower = if (
      model$k >= 3 && length(model$lower.predict) > 0
    ) {
      100 * plogis(first_or_na(model$lower.predict))
    } else {
      NA_real_
    },
    pi_upper = if (
      model$k >= 3 && length(model$upper.predict) > 0
    ) {
      100 * plogis(first_or_na(model$upper.predict))
    } else {
      NA_real_
    }
  )
}

collect_proportion_results <- function(model, analysis_name, group_name) {
  result <- extract_proportion_result(model)
  append_summary_row(
    analysis_name,
    "Proportion (%)",
    group_name,
    model$k,
    result$estimate,
    result$lower,
    result$upper,
    result$Q,
    result$df_Q,
    result$p_Q,
    result$I2,
    result$pi_lower,
    result$pi_upper
  )
}
