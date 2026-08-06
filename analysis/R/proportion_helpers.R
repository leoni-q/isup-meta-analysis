# Proportion model fitting, reporting, sensitivity, and plotting.

# Random-effects pooled proportion via GLMM (logit); k=1 groups reported exactly
fit_proportion <- function(
  source_data,
  event_column,
  total_column,
  title,
  include_subgroups = TRUE
) {
  has_event_and_total <- (
    !is.na(source_data[[event_column]]) &
      !is.na(source_data[[total_column]])
  )
  analysis_data <- source_data[has_event_and_total, ]

  model_arguments <- list(
    event = analysis_data[[event_column]],
    n = analysis_data[[total_column]],
    studlab = paste0(
      analysis_data$trial,
      " ",
      analysis_data$year,
      " (",
      analysis_data$arm,
      ")"
    ),
    sm = "PLOGIT",
    method = "GLMM",
    method.tau = "ML",
    random = TRUE,
    common = FALSE,
    title = title
  )

  if (include_subgroups) {
    model_arguments$subgroup <- droplevels(analysis_data$strategy)
    model_arguments$subgroup.name <- "Screening strategy"
  }

  suppressWarnings(do.call(metaprop, model_arguments))
}

# GLMM estimates are logits; plogis() converts them back to proportions.
report_proportion <- function(model, title) {
  log_section(title)
  log_line(
    sprintf(
      paste0(
        "k = %d arms | model: random effects GLMM (logit-normal), ",
        "back-transformed"
      ),
      model$k
    )
  )
  log_line(
    sprintf(
      "  Overall pooled proportion = %.4f%% (95%% CI %.4f%% to %.4f%%)",
      100 * plogis(model$TE.random),
      100 * plogis(model$lower.random),
      100 * plogis(model$upper.random)
    )
  )
  log_line(format_heterogeneity(model))

  if (!is.null(model$subgroup)) {
    log_line("  By screening strategy (random effects):")

    for (subgroup_index in seq_along(model$subgroup.levels)) {
      if (model$k.w[subgroup_index] == 1) {
        arm_index <- which(
          as.character(model$subgroup) == model$subgroup.levels[subgroup_index]
        )
        exact_confidence_interval <- binom.test(
          model$event[arm_index],
          model$n[arm_index]
        )$conf.int

        log_line(
          sprintf(
            paste0(
              "    %-28s k=1  %.4f%% ",
              "(exact 95%% CI %.4f%% to %.4f%%)  [%s]"
            ),
            model$subgroup.levels[subgroup_index],
            100 * model$event[arm_index] / model$n[arm_index],
            100 * exact_confidence_interval[1],
            100 * exact_confidence_interval[2],
            model$studlab[arm_index]
          )
        )
      } else {
        subgroup_q_p_value <- first_or_na(model$pval.Q.w[subgroup_index])
        subgroup_q_p_value_text <- if (is.na(subgroup_q_p_value)) {
          "NA"
        } else {
          format.pval(subgroup_q_p_value, digits = 3, eps = 1e-4)
        }

        log_line(
          sprintf(
            paste0(
              "    %-28s k=%d  %.4f%% ",
              "(95%% CI %.4f%% to %.4f%%)  Q p = %s"
            ),
            model$subgroup.levels[subgroup_index],
            model$k.w[subgroup_index],
            100 * plogis(model$TE.random.w[subgroup_index]),
            100 * plogis(model$lower.random.w[subgroup_index]),
            100 * plogis(model$upper.random.w[subgroup_index]),
            subgroup_q_p_value_text
          )
        )
      }
    }
  }

  log_line("  Per-arm observed proportions:")
  for (arm_index in seq_len(model$k.all)) {
    log_line(
      sprintf(
        "    %-26s %-28s %7d/%-7d = %.4f%%",
        model$studlab[arm_index],
        as.character(model$subgroup)[arm_index],
        model$event[arm_index],
        model$n[arm_index],
        100 * model$event[arm_index] / model$n[arm_index]
      )
    )
  }

  invisible(model)
}

report_leave_one_out_proportion <- function(model, title) {
  if (model$k < 3) {
    return(invisible(NULL))
  }

  leave_one_out_model <- suppressWarnings(
    metainf(model, pooled = "random")
  )
  log_line(sprintf("  Leave-one-out (random effects, %s):", title))

  for (arm_index in seq_len(model$k)) {
    log_line(
      sprintf(
        "    omitting %-26s %.4f%% (%.4f%% to %.4f%%)",
        leave_one_out_model$studlab[arm_index],
        100 * plogis(leave_one_out_model$TE[arm_index]),
        100 * plogis(leave_one_out_model$lower[arm_index]),
        100 * plogis(leave_one_out_model$upper[arm_index])
      )
    )
  }

  invisible(leave_one_out_model)
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
      leftlabs = c("Trial (arm)", "Events", "N", "%", "95% CI"),
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
