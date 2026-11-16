#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

Sys.setenv(
  REAL_DATA_DEFINITIONS_ONLY = "true",
  REAL_DATA_DESIGN = "adme",
  REAL_DATA_TARGET_FIT_CAP = "40"
)
source("02_run_pilot.R")

task_names <- c(
  "adme_caco2_related",
  "adme_caco2_heterogeneous"
)
regime_labels <- c(
  adme_caco2_related = "Same-assay source",
  adme_caco2_heterogeneous = "Different-assay source"
)

fit_interval_example <- function(task_name) {
  task_index <- which(tasks$task == task_name)
  stopifnot(length(task_index) == 1L)
  spec <- tasks[task_index]
  repetition <- 1L
  # Use the same prespecified target split in both source regimes so that the
  # two panels differ only in the auxiliary source supplied to the procedure.
  seed <- 20262115L
  data <- prepare_replication(spec, seed, p_keep = 300L)
  cfg <- make_real_config(data, seed)
  source_initials <- fit_source_initials(data$sources, cfg, seed + 10L)
  fitted <- fit_proposed(
    data$target_fit,
    data$target_cal,
    data$sources,
    cfg,
    source_initials,
    seed + 20L
  )

  prediction <- predict_linear_model(fitted$model, data$target_test$x)
  result <- data.table(
    task = task_name,
    source_regime = regime_labels[[task_name]],
    repetition = repetition,
    observed = data$target_test$y,
    prediction = prediction,
    lower = prediction - fitted$q,
    upper = prediction + fitted$q,
    selected = fitted$selected
  )
  result[, covered := observed >= lower & observed <= upper]
  setorder(result, prediction)
  result[, compound_rank := seq_len(.N)]
  result[, panel_label := source_regime]
  result
}

intervals <- rbindlist(lapply(task_names, fit_interval_example))
intervals[, panel_label := factor(
  panel_label,
  levels = unique(panel_label)
)]

figure_dir <- file.path("results", "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(intervals, file.path(
  "results", "adme_caco2_representative_intervals.csv"
))

p <- ggplot(intervals, aes(compound_rank)) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    fill = "#9ecae1",
    alpha = 0.55
  ) +
  geom_line(aes(y = prediction, colour = "Prediction"), linewidth = 0.65) +
  geom_line(aes(y = observed, colour = "Observed response"), linewidth = 0.45) +
  facet_wrap(~panel_label, ncol = 1) +
  scale_colour_manual(
    values = c("Prediction" = "#08519c", "Observed response" = "#cb181d"),
    name = NULL
  ) +
  labs(
    x = "Test compounds ordered by predicted response",
    y = "Standardized Caco-2 response"
  ) +
  theme_bw(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

ggsave(
  file.path(figure_dir, "adme_caco2_prediction_intervals.pdf"),
  p,
  width = 7.4,
  height = 5.2
)

print(intervals[, .(
  n_test = .N,
  coverage = mean(covered),
  mean_length = mean(upper - lower),
  selected = selected[1L]
), by = source_regime])
