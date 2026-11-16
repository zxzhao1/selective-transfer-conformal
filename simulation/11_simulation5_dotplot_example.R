args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1L) args[1L] else file.path(
  "simulation", "results", "simulation_5_full_pilotfixed_final200.csv"
)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
})

results <- fread(input_file)
method_order <- c(
  "Proposed", "F-AdaTrans", "Target-only", "Complete pooling",
  "TransLasso", "Ora-F-Ada"
)
scenario_order <- c("gaussian", "variance_shift", "heavy_tail")
scenario_labels <- c("Gaussian", "Variance shift", "Heavy-tailed\nsource errors")

summary <- melt(
  results,
  id.vars = c("scenario_label", "method"),
  measure.vars = c("coverage", "interval_length", "prediction_mse"),
  variable.name = "metric",
  value.name = "value"
)[, .(
  mean = mean(value),
  se = sd(value) / sqrt(.N)
), by = .(scenario_label, method, metric)]

summary[, method := factor(method, levels = method_order)]
summary[, scenario := factor(
  scenario_label, levels = scenario_order, labels = scenario_labels
)]

selection <- results[method == "Proposed", .(
  mean = mean(selected_transfer),
  se = sd(selected_transfer) / sqrt(.N)
), by = scenario_label]
selection[, scenario := factor(
  scenario_label, levels = scenario_order, labels = scenario_labels
)]

palette <- c(
  "Proposed" = "#0072B2", "F-AdaTrans" = "#D55E00",
  "Target-only" = "#009E73", "Complete pooling" = "#CC79A7",
  "TransLasso" = "#E69F00", "Ora-F-Ada" = "#000000"
)
shapes <- c(
  "Proposed" = 16, "F-AdaTrans" = 17, "Target-only" = 15,
  "Complete pooling" = 18, "TransLasso" = 8, "Ora-F-Ada" = 4
)

base_theme <- theme_classic(base_size = 9) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key.width = grid::unit(1.2, "lines"),
    axis.text.x = element_text(size = 8),
    plot.margin = margin(5, 8, 5, 5)
  )

make_metric_panel <- function(metric_name, y_label, reference = NULL) {
  panel <- ggplot(
    summary[metric == metric_name],
    aes(x = scenario, y = mean, color = method, shape = method)
  )
  if (!is.null(reference)) {
    panel <- panel + geom_hline(
      yintercept = reference, linetype = "dashed", color = "grey45",
      linewidth = 0.45
    )
  }
  panel +
    geom_errorbar(
      aes(ymin = mean - se, ymax = mean + se),
      position = position_dodge(width = 0.68), width = 0.08,
      linewidth = 0.4
    ) +
    geom_point(position = position_dodge(width = 0.68), size = 2.2) +
    scale_color_manual(values = palette, drop = FALSE) +
    scale_shape_manual(values = shapes, drop = FALSE) +
    labs(x = NULL, y = y_label) +
    base_theme
}

coverage_panel <- make_metric_panel(
  "coverage", "Empirical coverage", reference = 0.9
) + scale_y_continuous(limits = c(0.885, 0.915), breaks = c(0.89, 0.90, 0.91))
length_panel <- make_metric_panel(
  "interval_length", "Average interval length"
)
mse_panel <- make_metric_panel(
  "prediction_mse", "Mean squared prediction error"
)
selection_panel <- ggplot(selection, aes(x = scenario, y = mean)) +
  geom_errorbar(
    aes(ymin = pmax(mean - se, 0), ymax = pmin(mean + se, 1)),
    width = 0.08, linewidth = 0.45, color = palette[["Proposed"]]
  ) +
  geom_point(size = 2.4, shape = shapes[["Proposed"]],
             color = palette[["Proposed"]]) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = NULL, y = "Transfer-selection rate") +
  base_theme + theme(legend.position = "none")

legend <- get_legend(coverage_panel)
panels <- plot_grid(
  coverage_panel + theme(legend.position = "none"),
  length_panel + theme(legend.position = "none"),
  mse_panel + theme(legend.position = "none"),
  selection_panel,
  labels = c("A", "B", "C", "D"), label_size = 11,
  ncol = 2, align = "hv"
)
figure <- plot_grid(panels, legend, ncol = 1, rel_heights = c(1, 0.12))

output_base <- file.path(
  dirname(input_file), "figures",
  "simulation_5_full_pilotfixed_final200_dotplot_example"
)
ggsave(paste0(output_base, ".pdf"), figure, width = 7.2, height = 6.0,
       bg = "white")
ggsave(paste0(output_base, ".png"), figure, width = 7.2, height = 6.0,
       dpi = 320, bg = "white")
message("Saved example figure to ", output_base, ".{pdf,png}")

# Alternative layout: methods are labeled directly, and the three scenarios
# are separated within each method row. Selection rates remain in the table.
scenario_palette <- c(
  "Gaussian" = "#0072B2", "Variance shift" = "#D55E00",
  "Heavy-tailed\nsource errors" = "#009E73"
)
scenario_shapes <- c(
  "Gaussian" = 16, "Variance shift" = 17,
  "Heavy-tailed\nsource errors" = 15
)
forest_data <- copy(summary)
forest_data[, method_row := match(as.character(method), rev(method_order))]
forest_data[, scenario_offset := c(-0.22, 0, 0.22)[as.integer(scenario)]]
forest_data[, y_position := method_row + scenario_offset]

forest_theme <- theme_classic(base_size = 9) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(face = "bold", size = 10),
    plot.margin = margin(5, 8, 5, 5)
  )

make_forest_panel <- function(metric_name, x_label, panel_label,
                              x_limits, reference = NULL) {
  panel_data <- forest_data[metric == metric_name]
  panel <- ggplot(
    panel_data,
    aes(color = scenario, shape = scenario)
  )
  if (!is.null(reference)) {
    panel <- panel + geom_vline(
      xintercept = reference, linetype = "dashed", color = "grey45",
      linewidth = 0.45
    )
  }
  panel +
    geom_errorbar(
      aes(
        xmin = mean - se, xmax = mean + se,
        y = y_position
      ),
      width = 0.07, linewidth = 0.4, orientation = "y"
    ) +
    geom_point(aes(x = mean, y = y_position), size = 2.2) +
    scale_y_continuous(
      breaks = seq_along(method_order), labels = rev(method_order),
      limits = c(0.55, length(method_order) + 0.45)
    ) +
    scale_x_continuous(limits = x_limits) +
    scale_color_manual(values = scenario_palette) +
    scale_shape_manual(values = scenario_shapes) +
    labs(x = x_label, title = panel_label) +
    forest_theme
}

forest_coverage <- make_forest_panel(
  "coverage", "Empirical coverage", "A", c(0.885, 0.915), 0.9
)
forest_length <- make_forest_panel(
  "interval_length", "Average interval length", "B", c(3.0, 7.2)
)
forest_mse <- make_forest_panel(
  "prediction_mse", "Mean squared prediction error", "C", c(0.9, 4.6)
)

forest_legend <- get_legend(forest_coverage)
forest_panels <- plot_grid(
  forest_coverage + theme(legend.position = "none"),
  forest_length + theme(legend.position = "none"),
  forest_mse + theme(legend.position = "none"),
  nrow = 1, align = "hv"
)
forest_figure <- plot_grid(
  forest_panels, forest_legend, ncol = 1, rel_heights = c(1, 0.13)
)
forest_output <- file.path(
  dirname(input_file), "figures",
  "simulation_5_full_pilotfixed_final200_forest_example"
)
ggsave(paste0(forest_output, ".pdf"), forest_figure,
       width = 8.2, height = 3.8, bg = "white")
ggsave(paste0(forest_output, ".png"), forest_figure,
       width = 8.2, height = 3.8, dpi = 320, bg = "white")
message("Saved forest example to ", forest_output, ".{pdf,png}")
