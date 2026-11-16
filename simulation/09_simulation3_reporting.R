#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)
script_path <- if (length(file_arg) == 1L) {
  normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
} else {
  normalizePath("simulation/09_simulation3_reporting.R", mustWork = TRUE)
}
simulation_dir <- dirname(script_path)
results_dir <- file.path(simulation_dir, "results")
summary_path <- file.path(
  results_dir,
  "tables",
  "simulation_3_full_fixed_final200_summary_long.csv"
)
figure_dir <- file.path(results_dir, "figures")

if (!file.exists(summary_path)) {
  stop("Simulation 3 summary not found: ", summary_path, call. = FALSE)
}

summary_data <- fread(summary_path)
plot_data <- summary_data[
  metric %chin% c("interval_length", "prediction_mse", "selection_rate")
]

plot_data[, experiment := fcase(
  grepl("^sim3_nsource_", scenario), "Source sample size",
  grepl("^sim3_K_", scenario), "Number of sources",
  grepl("^sim3_nfit_", scenario), "Target fitting sample size",
  default = NA_character_
)]
plot_data[, factor_value := as.numeric(sub(".*_", "", scenario))]
plot_data[, factor_axis := factor(
  factor_value,
  levels = c(2, 4, 8, 30, 50, 100, 250, 500)
)]
plot_data[, metric_label := factor(
  metric,
  levels = c("interval_length", "prediction_mse", "selection_rate"),
  labels = c(
    "Average interval length",
    "Prediction MSE",
    "Transfer-selection rate"
  )
)]
plot_data[, experiment := factor(
  experiment,
  levels = c(
    "Source sample size",
    "Number of sources",
    "Target fitting sample size"
  )
)]
plot_data[, method := factor(
  method,
  levels = c("Proposed", "F-AdaTrans", "Target-only")
)]

if (anyNA(plot_data$experiment) || anyNA(plot_data$factor_value)) {
  stop("Simulation 3 scenario names do not match the expected design.",
       call. = FALSE)
}

palette <- c(
  "Proposed" = "#0072B2",
  "F-AdaTrans" = "#D55E00",
  "Target-only" = "#009E73"
)

result_figure <- ggplot(
  plot_data,
  aes(
    x = factor_axis,
    y = mean,
    color = method,
    shape = method,
    group = method
  )
) +
  geom_line(linewidth = 0.55, na.rm = TRUE) +
  geom_point(size = 2.2, na.rm = TRUE) +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    width = 0,
    linewidth = 0.45,
    na.rm = TRUE
  ) +
  facet_grid(
    rows = vars(metric_label),
    cols = vars(experiment),
    scales = "free",
    switch = "y"
  ) +
  scale_color_manual(values = palette, drop = FALSE) +
  scale_x_discrete(drop = TRUE) +
  labs(x = NULL, y = NULL, color = "Method", shape = "Method") +
  theme_bw(base_size = 10.5) +
  theme(
    legend.position = "bottom",
    legend.box.margin = margin(t = -4),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = "grey65"),
    strip.placement = "outside",
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    panel.spacing = grid::unit(0.8, "lines")
  )

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
figure_base <- file.path(
  figure_dir,
  "simulation_3_full_fixed_final200_faceted_results"
)
ggsave(
  paste0(figure_base, ".pdf"),
  result_figure,
  width = 9.2,
  height = 7.2,
  device = grDevices::pdf
)
ggsave(
  paste0(figure_base, ".png"),
  result_figure,
  width = 9.2,
  height = 7.2,
  dpi = 320
)

message("Wrote Simulation 3 figure to ", figure_base, ".{pdf,png}")
