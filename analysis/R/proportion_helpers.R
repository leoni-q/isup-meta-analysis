# Proportion model fitting, reporting, sensitivity, and plotting.

# Random-effects pooled proportion via GLMM (logit); k=1 groups reported exactly
fit_proportion <- function(
  source_data,
  event_column,
  total_column,
  title,
  include_subgroups = FALSE
) {
  has_event_and_total <- (
    !is.na(source_data[[event_column]]) &
      !is.na(source_data[[total_column]])
  )
  analysis_data <- source_data[has_event_and_total, ]
  stopifnot(
    nrow(analysis_data) > 0,
    all(analysis_data[[event_column]] <= analysis_data[[total_column]])
  )

  source_arm_labels <- if ("source_arms" %in% names(analysis_data)) {
    analysis_data$source_arms
  } else {
    analysis_data$arm
  }

  model_arguments <- list(
    event = analysis_data[[event_column]],
    n = analysis_data[[total_column]],
    studlab = paste0(
      analysis_data$trial,
      " ",
      analysis_data$year,
      " (",
      source_arm_labels,
      ")"
    ),
    sm = "PLOGIT",
    method = "GLMM",
    method.tau = "ML",
    random = TRUE,
    common = FALSE,
    prediction = TRUE,
    title = title
  )

  if (include_subgroups) {
    model_arguments$subgroup <- droplevels(analysis_data$strategy)
    model_arguments$subgroup.name <- "Screening strategy"
  }

  suppressWarnings(do.call(metaprop, model_arguments))
}

# GLMM estimates are logits; plogis() converts them back to proportions.
report_proportion <- function(model, title, group_name) {
  log_section(title)
  log_line(
    sprintf(
      paste0(
        "k = %d studies | model: random effects GLMM (logit-normal), ",
        "back-transformed"
      ),
      model$k
    )
  )
  result <- extract_proportion_result(model)
  interval_label <- if (model$k == 1) "exact 95% CI" else "95% CI"
  log_line(
    sprintf(
      "  %s: %.4f%% (%s %.4f%% to %.4f%%)",
      group_name,
      result$estimate,
      interval_label,
      result$lower,
      result$upper
    )
  )
  if (model$k > 1) {
    log_line(format_heterogeneity(model))
  }
  if (!is.na(result$pi_lower)) {
    log_line(
      sprintf(
        "  95%% prediction interval: %.4f%% to %.4f%%",
        result$pi_lower,
        result$pi_upper
      )
    )
  }

  log_line("  Per-study observed proportions:")
  for (study_index in seq_len(model$k.all)) {
    log_line(
      sprintf(
        "    %-28s %7d/%-7d = %.4f%%",
        model$studlab[study_index],
        model$event[study_index],
        model$n[study_index],
        100 * model$event[study_index] / model$n[study_index]
      )
    )
  }

  invisible(model)
}

report_leave_one_out_proportion <- function(
  source_data,
  event_column,
  total_column,
  analysis_name,
  group_name
) {
  has_event_and_total <- (
    !is.na(source_data[[event_column]]) &
      !is.na(source_data[[total_column]])
  )
  analysis_data <- source_data[has_event_and_total, ]

  if (nrow(analysis_data) < 2) {
    log_line(
      sprintf(
        "  Leave-one-study-out skipped for %s / %s (k = %d < 2).",
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

  for (study_index in seq_len(nrow(analysis_data))) {
    omitted_study <- analysis_data$trial[study_index]
    remaining_data <- analysis_data[-study_index, , drop = FALSE]
    sensitivity_model <- suppressWarnings(
      fit_proportion(
        source_data = remaining_data,
        event_column = event_column,
        total_column = total_column,
        title = paste0(analysis_name, " omitting ", omitted_study),
        include_subgroups = FALSE
      )
    )
    leave_one_out_models[[study_index]] <- sensitivity_model
    result <- extract_proportion_result(sensitivity_model)

    log_line(
      sprintf(
        "    omitting %-22s %.4f%% (%.4f%% to %.4f%%)",
        omitted_study,
        result$estimate,
        result$lower,
        result$upper
      )
    )
    append_leave_one_out_row(
      analysis = analysis_name,
      group = group_name,
      omitted_study = omitted_study,
      remaining_studies = sensitivity_model$k,
      scale = "Proportion (%)",
      estimate = result$estimate,
      confidence_interval_lower = result$lower,
      confidence_interval_upper = result$upper,
      i_squared_percent = result$I2
    )
  }

  invisible(leave_one_out_models)
}

save_proportion_forest_plot <- function(
  model,
  file_stem,
  x_axis_label,
  summary_label = "Pooled proportion (%)"
) {
  subgroup_count <- if (is.null(model$subgroup)) {
    0
  } else {
    length(model$subgroup.levels)
  }
  plot_height <- 0.8 + 0.24 * (model$k.all + 3 * subgroup_count + 5)

  for (file_format in c("pdf", "png")) {
    output_path <- file.path(
      OUTPUT_DIR,
      paste0(file_stem, ".", file_format)
    )

    if (file_format == "pdf") {
      pdf(output_path, width = 11, height = plot_height)
    } else {
      png(
        output_path,
        width = 11,
        height = plot_height,
        units = "in",
        res = 200
      )
    }

    forest(
      model,
      layout = "meta",
      common = FALSE,
      random = TRUE,
      xlab = x_axis_label,
      pscale = 100,
      smlab = summary_label,
      leftcols = c("studlab", "event", "n", "effect", "ci"),
      leftlabs = c("Study (source arms)", "Events", "N", "%", "95% CI"),
      rightcols = FALSE,
      print.subgroup.name = FALSE,
      col.subgroup = "black",
      print.I2 = TRUE,
      print.Q = TRUE,
      print.pval.Q = TRUE,
      digits = 2,
      col.diamond = "#b2182b",
      col.diamond.lines = "black",
      col.square = "#d6604d",
      col.square.lines = "black"
    )
    grid::grid.text(
      model$title,
      x = 0.5,
      y = 0.99,
      just = "top",
      gp = grid::gpar(cex = 1.0, fontface = "bold")
    )
    dev.off()
  }
}
