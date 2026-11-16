args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Provide the Simulation 4 result CSV.")

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

results <- fread(args[1L])
metrics <- c("interval_length", "prediction_mse")
summary <- melt(
  results,
  id.vars = c("scenario_number", "scenario_label", "method"),
  measure.vars = metrics,
  variable.name = "metric",
  value.name = "value"
)[, .(
  mean = mean(value),
  se = sd(value) / sqrt(.N)
), by = .(scenario_number, scenario_label, method, metric)]

method_order <- c(
  "Proposed", "F-AdaTrans", "Target-only", "Complete pooling",
  "TransLasso", "Ora-F-Ada"
)
summary[, method := factor(method, levels = method_order)]
scenario_levels <- sort(unique(summary$scenario_number))
scenario_labels <- if (length(scenario_levels) == 3L) {
  c("Gaussian", "Variance shift", "Heavy-tailed\nsource errors")
} else {
  c(
    "Complete\ntransferability",
    "Sparse feature-wise\nheterogeneity",
    "Heterogeneous source\ntransferability",
    "Dense feature-wise\nheterogeneity"
  )
}
summary[, scenario := factor(
  scenario_number, levels = scenario_levels, labels = scenario_labels
)]
palette <- c(
  "Proposed" = "#0072B2", "F-AdaTrans" = "#D55E00",
  "Target-only" = "#009E73", "Complete pooling" = "#CC79A7",
  "TransLasso" = "#E69F00", "Ora-F-Ada" = "#000000"
)

figure_dir <- file.path(dirname(args[1L]), "figures")
base_name <- tools::file_path_sans_ext(basename(args[1L]))
for (metric_name in metrics) {
  plot_data <- summary[metric == metric_name]
  figure <- ggplot(
    plot_data,
    aes(x = scenario, y = mean, color = method, shape = method, group = method)
  ) +
    geom_line(linewidth = 0.55, na.rm = TRUE) +
    geom_point(size = 2.2, na.rm = TRUE) +
    geom_errorbar(
      aes(ymin = pmax(mean - se, .Machine$double.eps), ymax = mean + se),
      width = 0.1, linewidth = 0.4, na.rm = TRUE
    ) +
    scale_y_log10() +
    scale_color_manual(values = palette, drop = FALSE) +
    labs(
      x = NULL,
      y = if (metric_name == "prediction_mse") {
        "Mean squared prediction error (log scale)"
      } else {
        "Average interval length (log scale)"
      },
      color = "Method", shape = "Method"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 9)
    )
  output_base <- file.path(
    figure_dir, paste0(base_name, "_", metric_name, "_log")
  )
  ggsave(paste0(output_base, ".pdf"), figure, width = 7.2, height = 4.8)
  ggsave(paste0(output_base, ".png"), figure, width = 7.2, height = 4.8,
         dpi = 320)
}

selection <- results[method == "Proposed", .(
  mean = mean(selected_transfer),
  se = sd(selected_transfer) / sqrt(.N)
), by = scenario_number]
selection[, scenario := factor(
  scenario_number, levels = scenario_levels, labels = scenario_labels
)]
selection_figure <- ggplot(
  selection, aes(x = scenario, y = mean, group = 1)
) +
  geom_line(color = "#0072B2", linewidth = 0.6) +
  geom_point(color = "#0072B2", size = 2.4) +
  geom_errorbar(
    aes(ymin = pmax(mean - se, 0), ymax = pmin(mean + se, 1)),
    width = 0.1, color = "#0072B2"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = "Transfer-selection rate") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(size = 9))
output_base <- file.path(
  figure_dir, paste0(base_name, "_selection_rate_compact")
)
ggsave(paste0(output_base, ".pdf"), selection_figure, width = 7.2, height = 4.8)
ggsave(paste0(output_base, ".png"), selection_figure, width = 7.2, height = 4.8,
       dpi = 320)
