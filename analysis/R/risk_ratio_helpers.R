# Risk-ratio model fitting, reporting, sensitivity, and plotting.

# Random-effects RR: inverse variance + DL tau + Hartung-Knapp CI, incr 0.5
# The meta package calls intervention fields "experimental" (.e) and
# comparator fields "control" (.c), hence names such as event.e and event.c.
complete_rr_data <- function(
  source_data,
  intervention_event_column,
  intervention_total_column,
  comparator_event_column,
  comparator_total_column
) {
  is_complete <- (
    !is.na(source_data[[intervention_event_column]]) &
      !is.na(source_data[[intervention_total_column]]) &
      !is.na(source_data[[comparator_event_column]]) &
      !is.na(source_data[[comparator_total_column]])
  )
  source_data[is_complete, , drop = FALSE]
}

fit_risk_ratio <- function(
  source_data,
  intervention_event_column,
  intervention_total_column,
  comparator_event_column,
  comparator_total_column,
  title,
  include_subgroups = TRUE
) {
  analysis_data <- complete_rr_data(
    source_data,
    intervention_event_column,
    intervention_total_column,
    comparator_event_column,
    comparator_total_column
  )
  stopifnot(
    nrow(analysis_data) > 0,
    all(analysis_data[[intervention_event_column]] <=
      analysis_data[[intervention_total_column]]),
    all(analysis_data[[comparator_event_column]] <=
      analysis_data[[comparator_total_column]])
  )

  model_arguments <- list(
    event.e = analysis_data[[intervention_event_column]],
    n.e = analysis_data[[intervention_total_column]],
    event.c = analysis_data[[comparator_event_column]],
    n.c = analysis_data[[comparator_total_column]],
    studlab = analysis_data$label,
    sm = "RR",
    method = "Inverse",
    method.tau = "DL",
    method.random.ci = "HK",
    incr = 0.5,
    method.incr = "only0",
    RR.Cochrane = TRUE,
    random = TRUE,
    common = FALSE,
    prediction = TRUE,
    title = title
  )

  if (
    include_subgroups &&
      length(unique(analysis_data$family_label)) > 1
  ) {
    model_arguments$subgroup <- analysis_data$family_label
    model_arguments$subgroup.name <- "Contrast family"
    model_arguments$test.subgroup <- TRUE
  }

  do.call(metabin, model_arguments)
}

# `meta` stores risk ratios on the log scale; exp() converts them back to RRs.
report_risk_ratio <- function(model, title) {
  log_section(title)
  # With a single contrast there is nothing to pool: `meta` reports the study's
  # own log-RR with a normal-approximation interval, so neither the DL tau
  # estimator nor the Hartung-Knapp adjustment is in effect.
  log_line(
    if (model$k == 1) {
      paste0(
        "k = 1 trial | single-study contrast: no pooling, no tau estimated; ",
        "CI is the normal-approximation (Wald) interval on the log scale"
      )
    } else {
      sprintf(
        paste0(
          "k = %d trials | model: random effects, inverse variance, ",
          "DL tau, Hartung-Knapp CI"
        ),
        model$k
      )
    }
  )
  log_line(
    sprintf(
      "  Pooled RR = %.3f (95%% CI %.3f to %.3f), p = %s",
      exp(model$TE.random),
      exp(model$lower.random),
      exp(model$upper.random),
      format.pval(model$pval.random, digits = 4, eps = 1e-4)
    )
  )
  if (model$k > 1) {
    log_line(format_heterogeneity(model))
  }
  if (
    model$k >= 3 &&
      length(model$lower.predict) > 0 &&
      !is.na(first_or_na(model$lower.predict))
  ) {
    log_line(
      sprintf(
        "  95%% prediction interval: %.3f to %.3f",
        exp(first_or_na(model$lower.predict)),
        exp(first_or_na(model$upper.predict))
      )
    )
  }

  if (!is.null(model$subgroup)) {
    log_line("  Subgroup estimates (random effects):")

    for (subgroup_index in seq_along(model$subgroup.levels)) {
      log_line(
        sprintf(
          "    %-55s k=%d  RR = %.3f (%.3f to %.3f)",
          model$subgroup.levels[subgroup_index],
          model$k.w[subgroup_index],
          exp(model$TE.random.w[subgroup_index]),
          exp(model$lower.random.w[subgroup_index]),
          exp(model$upper.random.w[subgroup_index])
        )
      )
    }

    if (!is.na(model$pval.Q.b.random)) {
      log_line(
        sprintf(
          "    Test for subgroup differences: Q = %.2f, df = %d, p = %s",
          model$Q.b.random,
          model$df.Q.b,
          format.pval(model$pval.Q.b.random, digits = 4, eps = 1e-4)
        )
      )
    }
  }

  log_line("  Per-trial risk ratios:")
  for (trial_index in seq_len(model$k.all)) {
    log_line(
      sprintf(
        "    %-22s %6d/%-7d vs %6d/%-7d   RR = %.3f (%.3f to %.3f)",
        model$studlab[trial_index],
        model$event.e[trial_index],
        model$n.e[trial_index],
        model$event.c[trial_index],
        model$n.c[trial_index],
        exp(model$TE[trial_index]),
        exp(model$lower[trial_index]),
        exp(model$upper[trial_index])
      )
    )
  }

  invisible(model)
}

report_leave_one_out_risk_ratio <- function(
  source_data,
  intervention_event_column,
  intervention_total_column,
  comparator_event_column,
  comparator_total_column,
  analysis_name,
  group_name
) {
  analysis_data <- complete_rr_data(
    source_data,
    intervention_event_column,
    intervention_total_column,
    comparator_event_column,
    comparator_total_column
  )

  # Omitting one of two studies leaves a single trial, whose own estimate is not
  # a pooled result; publishing it as one would misrepresent the sensitivity.
  if (nrow(analysis_data) < 3) {
    log_line(
      sprintf(
        paste0(
          "  Leave-one-study-out skipped for %s / %s (k = %d < 3; omitting a ",
          "study would leave no pooled estimate)."
        ),
        analysis_name,
        group_name,
        nrow(analysis_data)
      )
    )
    return(invisible(NULL))
  }

  log_line(
    sprintf(
      "  Leave-one-study-out (%s; %s):",
      analysis_name,
      group_name
    )
  )

  leave_one_out_models <- vector("list", nrow(analysis_data))

  for (trial_index in seq_len(nrow(analysis_data))) {
    omitted_trial <- analysis_data$trial[trial_index]
    remaining_data <- analysis_data[-trial_index, , drop = FALSE]
    sensitivity_model <- suppressWarnings(
      fit_risk_ratio(
        source_data = remaining_data,
        intervention_event_column = intervention_event_column,
        intervention_total_column = intervention_total_column,
        comparator_event_column = comparator_event_column,
        comparator_total_column = comparator_total_column,
        title = paste0(analysis_name, " omitting ", omitted_trial),
        include_subgroups = FALSE
      )
    )
    leave_one_out_models[[trial_index]] <- sensitivity_model

    i_squared <- if (sensitivity_model$k > 1) {
      100 * first_or_na(sensitivity_model$I2)
    } else {
      NA_real_
    }
    i_squared_text <- if (is.na(i_squared)) {
      "NA"
    } else {
      sprintf("%.1f%%", i_squared)
    }

    log_line(
      sprintf(
        "    omitting %-22s RR = %.3f (%.3f to %.3f)  I2 = %s",
        omitted_trial,
        exp(sensitivity_model$TE.random),
        exp(sensitivity_model$lower.random),
        exp(sensitivity_model$upper.random),
        i_squared_text
      )
    )
    append_leave_one_out_row(
      analysis = analysis_name,
      group = group_name,
      omitted_study = omitted_trial,
      remaining_studies = sensitivity_model$k,
      scale = "RR",
      estimate = exp(sensitivity_model$TE.random),
      confidence_interval_lower = exp(sensitivity_model$lower.random),
      confidence_interval_upper = exp(sensitivity_model$upper.random),
      i_squared_percent = i_squared
    )
  }

  invisible(leave_one_out_models)
}

# Preferred ticks are tuned for the pooled cohorts; single-trial contrasts such
# as PROBASE can sit far outside them. Widening the axis to cover every point
# estimate keeps squares and diamonds on the plot - wide confidence intervals
# are still drawn as arrows, which is the conventional forest-plot behaviour.
# Progressively wider tick sets. The first one that spans every point estimate
# is used, so the axis stays as tight as the data allows while never cutting an
# estimate off. Each set includes 1 so the null line is always drawn.
FOREST_TICK_SETS <- list(
  c(0.2, 0.35, 0.5, 0.7, 1, 1.5, 2, 3),
  c(0.2, 0.5, 1, 2, 5),
  c(0.1, 0.2, 0.5, 1, 2, 5, 10),
  c(0.5, 1, 2, 5, 10, 25),
  c(0.1, 0.5, 1, 2, 5, 10, 50),
  c(0.1, 0.5, 1, 5, 20, 100),
  c(0.1, 1, 10, 100, 500),
  c(0.01, 0.1, 1, 10, 100, 1000)
)

risk_ratio_x_ticks <- function(model, preferred_ticks) {
  point_estimates <- exp(c(model$TE, model$TE.random))
  point_estimates <- point_estimates[
    is.finite(point_estimates) & point_estimates > 0
  ]
  if (length(point_estimates) == 0) {
    return(preferred_ticks)
  }

  covers_estimates <- function(ticks) {
    min(point_estimates) >= min(ticks) && max(point_estimates) <= max(ticks)
  }
  if (covers_estimates(preferred_ticks)) {
    return(preferred_ticks)
  }
  for (candidate_ticks in FOREST_TICK_SETS) {
    if (covers_estimates(candidate_ticks)) {
      return(candidate_ticks)
    }
  }

  # Beyond the widest prepared set, fall back to a padded log-spaced axis.
  lower_limit <- min(c(point_estimates, 1)) / 2
  upper_limit <- max(c(point_estimates, 1)) * 2
  signif(exp(seq(log(lower_limit), log(upper_limit), length.out = 5)), 1)
}

save_risk_ratio_forest_plot <- function(
  model,
  file_stem,
  left_label,
  right_label,
  x_ticks = c(0.2, 0.5, 1, 2, 5, 10)
) {
  x_ticks <- risk_ratio_x_ticks(model, x_ticks)
  subgroup_count <- if (is.null(model$subgroup)) {
    0
  } else {
    length(model$subgroup.levels)
  }
  plot_height <- 0.8 + 0.28 * (model$k.all + 3 * subgroup_count + 5)

  for (file_format in c("pdf", "png")) {
    output_path <- file.path(
      OUTPUT_DIR,
      paste0(file_stem, ".", file_format)
    )

    if (file_format == "pdf") {
      pdf(output_path, width = 11.5, height = plot_height)
    } else {
      png(
        output_path,
        width = 11.5,
        height = plot_height,
        units = "in",
        res = 200
      )
    }

    forest(
      model,
      layout = "RevMan5",
      common = FALSE,
      random = TRUE,
      label.e = "Intervention",
      label.c = "Comparator",
      label.left = left_label,
      label.right = right_label,
      xlim = range(x_ticks),
      at = x_ticks,
      print.subgroup.name = FALSE,
      col.subgroup = "black",
      # Standard heterogeneity line disabled: it labels Q as Chi^2 and
      # reports I^2. Replaced by a Cochran's Q line without I^2.
      hetstat = FALSE,
      text.addline1 = format_forest_heterogeneity_label(model),
      fs.addline = 10,
      digits = 2,
      col.diamond = "#2166ac",
      col.diamond.lines = "black",
      col.square = "#4393c3",
      col.square.lines = "black"
    )
    grid::grid.text(
      model$title,
      x = 0.5,
      y = 0.99,
      just = "top",
      gp = grid::gpar(cex = 1.05, fontface = "bold")
    )
    dev.off()
  }
}
