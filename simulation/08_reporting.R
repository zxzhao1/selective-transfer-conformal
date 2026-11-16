monte_carlo_summary <- function(results) {
  metrics <- c(
    "coverage", "interval_length", "prediction_mse", "slope_error"
  )
  id_columns <- c(
    "simulation", "scenario", "scenario_label", "method"
  )
  missing <- setdiff(c(id_columns, metrics, "selected_transfer"), names(results))
  if (length(missing) > 0L) {
    stop("Results are missing columns: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }

  data <- data.table::as.data.table(results)
  scenario_order <- unique(data$scenario)
  label_order <- unique(data$scenario_label)
  long <- data.table::melt(
    data,
    id.vars = id_columns,
    measure.vars = metrics,
    variable.name = "metric",
    value.name = "value"
  )
  summary <- long[, .(
    mean = mean(value, na.rm = TRUE),
    se = if (.N > 1L) stats::sd(value, na.rm = TRUE) / sqrt(.N) else NA_real_,
    replications = .N
  ), by = c(id_columns, "metric")]

  selection <- data[!is.na(selected_transfer), .(
    mean = mean(selected_transfer),
    se = if (.N > 1L) stats::sd(selected_transfer) / sqrt(.N) else NA_real_,
    replications = .N
  ), by = id_columns]
  if (nrow(selection) > 0L) {
    selection[, metric := "selection_rate"]
    summary <- data.table::rbindlist(
      list(summary, selection), use.names = TRUE, fill = TRUE
    )
  }

  summary[, scenario := factor(scenario, levels = scenario_order)]
  summary[, scenario_label := factor(scenario_label, levels = label_order)]
  summary[]
}

format_mean_se <- function(mean, se, digits = 3L) {
  mean_text <- formatC(mean, format = "f", digits = digits)
  ifelse(
    is.na(se),
    paste0(mean_text, " (--)"),
    paste0(mean_text, " (", formatC(se, format = "f", digits = digits), ")")
  )
}

escape_latex <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("_", "\\\\_", x, fixed = TRUE)
  x <- gsub("%", "\\\\%", x, fixed = TRUE)
  x <- gsub("&", "\\\\&", x, fixed = TRUE)
  x
}

format_scenario_latex <- function(x) {
  if (grepl("^[hs]_delta=", x)) {
    symbol <- sub("_delta=.*$", "", x)
    value <- sub("^[hs]_delta=", "", x)
    return(paste0("$", symbol, "_{\\delta}=", value, "$"))
  }
  escape_latex(x)
}

write_latex_summary <- function(summary, path) {
  table_data <- data.table::copy(summary)
  table_data[, value := format_mean_se(mean, se)]
  wide <- data.table::dcast(
    table_data,
    simulation + scenario + scenario_label + method ~ metric,
    value.var = "value",
    fill = "--"
  )
  desired <- c(
    "coverage", "interval_length", "prediction_mse", "slope_error",
    "selection_rate"
  )
  for (column in setdiff(desired, names(wide))) {
    wide[, (column) := "--"]
  }

  lines <- c(
    "\\begin{table}[htbp]",
    "\\centering",
    "\\caption{Monte Carlo performance. Entries are means with Monte Carlo standard errors in parentheses.}",
    "\\label{tab:simulation-results}",
    "\\begin{tabular}{llccccc}",
    "\\hline",
    "Scenario & Method & Coverage & Length & Prediction MSE & Slope error & Selection rate \\\\",
    "\\hline"
  )
  for (row in seq_len(nrow(wide))) {
    lines <- c(lines, paste0(
      format_scenario_latex(as.character(wide$scenario_label[row])), " & ",
      escape_latex(wide$method[row]), " & ",
      wide$coverage[row], " & ",
      wide$interval_length[row], " & ",
      wide$prediction_mse[row], " & ",
      wide$slope_error[row], " & ",
      wide$selection_rate[row], " \\\\"
    ))
  }
  lines <- c(lines, "\\hline", "\\end{tabular}", "\\end{table}")
  writeLines(lines, path)
  invisible(wide)
}

metric_axis_label <- function(metric) {
  switch(
    metric,
    coverage = "Empirical coverage",
    interval_length = "Average interval length",
    prediction_mse = "Mean squared prediction error",
    slope_error = "Slope estimation error",
    selection_rate = "Transfer-selection rate",
    metric
  )
}

make_metric_plot <- function(summary, metric_name, alpha) {
  plot_data <- as.data.frame(summary)
  plot_data <- plot_data[
    as.character(plot_data$metric) == metric_name, , drop = FALSE
  ]
  original_labels <- as.character(plot_data$scenario_label)
  parse_delta_parameter <- all(grepl("^[hs]_delta=", original_labels))
  if (parse_delta_parameter) {
    delta_symbol <- sub("_delta=.*$", "", original_labels[1L])
    label_levels <- unique(original_labels)
    parsed_levels <- paste0(
      delta_symbol, "[delta] == ", sub("^[hs]_delta=", "", label_levels)
    )
    plot_data$scenario_axis <- factor(
      paste0(
        delta_symbol, "[delta] == ",
        sub("^[hs]_delta=", "", original_labels)
      ),
      levels = parsed_levels
    )
  } else {
    plot_data$scenario_axis <- plot_data$scenario_label
  }
  plot_data$se_plot <- ifelse(is.na(plot_data$se), 0, plot_data$se)
  palette <- c(
    "Proposed" = "#0072B2",
    "F-AdaTrans" = "#D55E00",
    "Target-only" = "#009E73",
    "Complete pooling" = "#CC79A7",
    "TransLasso" = "#E69F00",
    "Ora-F-Ada" = "#000000"
  )

  figure <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = scenario_axis,
      y = mean,
      color = method,
      group = method,
      shape = method
    )
  ) +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean - se_plot, ymax = mean + se_plot),
      width = 0.12,
      linewidth = 0.45
    ) +
    ggplot2::scale_color_manual(values = palette, drop = FALSE) +
    ggplot2::labs(
      x = NULL,
      y = metric_axis_label(metric_name),
      color = "Method",
      shape = "Method"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
      panel.grid.minor = ggplot2::element_blank()
    )

  if (length(unique(plot_data$scenario_label)) > 1L) {
    figure <- figure + ggplot2::geom_line(linewidth = 0.6)
  }

  if (parse_delta_parameter) {
    figure <- figure + ggplot2::scale_x_discrete(
      labels = scales::parse_format()
    )
  }

  if (metric_name == "coverage") {
    figure <- figure +
      ggplot2::geom_hline(
        yintercept = 1 - alpha,
        linetype = "dashed",
        color = "grey35"
      ) +
      ggplot2::coord_cartesian(ylim = c(0, 1))
  }
  if (metric_name == "selection_rate") {
    figure <- figure + ggplot2::coord_cartesian(ylim = c(0, 1))
  }
  figure
}

write_simulation_outputs <- function(results, output_stem, alpha = 0.1) {
  summary <- monte_carlo_summary(results)
  output_dir <- dirname(output_stem)
  table_dir <- file.path(output_dir, "tables")
  figure_dir <- file.path(output_dir, "figures")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

  base_name <- basename(output_stem)
  long_path <- file.path(table_dir, paste0(base_name, "_summary_long.csv"))
  wide_path <- file.path(table_dir, paste0(base_name, "_summary_wide.csv"))
  tex_path <- file.path(table_dir, paste0(base_name, "_summary.tex"))
  data.table::fwrite(summary, long_path)

  wide <- data.table::dcast(
    summary,
    simulation + scenario + scenario_label + method ~ metric,
    value.var = c("mean", "se", "replications")
  )
  data.table::fwrite(wide, wide_path)
  write_latex_summary(summary, tex_path)

  metrics <- unique(as.character(summary$metric))
  for (metric_name in metrics) {
    figure <- make_metric_plot(summary, metric_name, alpha)
    file_base <- file.path(figure_dir, paste0(base_name, "_", metric_name))
    ggplot2::ggsave(
      paste0(file_base, ".pdf"), figure, width = 7.2, height = 4.8,
      device = grDevices::pdf
    )
    ggplot2::ggsave(
      paste0(file_base, ".png"), figure, width = 7.2, height = 4.8,
      dpi = 320
    )
  }

  list(
    summary = summary,
    table_long = long_path,
    table_wide = wide_path,
    table_tex = tex_path,
    figure_dir = figure_dir
  )
}
