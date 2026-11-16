#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(gridExtra)
})

root <- normalizePath(getwd(), mustWork = TRUE)
result_dir <- file.path(root, "results")
figure_dir <- file.path(result_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

main <- fread(file.path(
  result_dir, "real_results_adme_nfit40_10rep_p300.csv"
))
translasso <- fread(file.path(
  result_dir, "real_results_translasso_adme_nfit40_10rep_p300.csv"
))
results <- rbind(main, translasso, fill = TRUE)
tasks <- fread(file.path(root, "data", "processed", "adme_tasks.csv"))
task_summary <- fread(file.path(
  result_dir, "adme_nfit40_10rep_task_summary.csv"
))

metric_columns <- c(
  "coverage", "interval_length", "prediction_mse", "prediction_mae"
)
within_task_sd <- results[, lapply(.SD, sd), by = .(method, task),
                          .SDcols = metric_columns]
stability_summary <- within_task_sd[, lapply(.SD, mean), by = method,
                                    .SDcols = metric_columns]
setnames(
  stability_summary,
  metric_columns,
  paste0(metric_columns, "_mean_within_task_sd")
)
fwrite(
  stability_summary,
  file.path(result_dir, "adme_nfit40_10rep_stability_summary.csv")
)

endpoint_labels <- c(
  caco2 = "Caco-2 permeability",
  vdss = "VDss",
  ppbr = "PPBR",
  half_life = "Half-life",
  clearance_hepatocyte = "Hepatocyte clearance",
  clearance_microsome = "Microsomal clearance",
  solubility = "Aqueous solubility",
  lipophilicity = "Lipophilicity",
  hydration = "Hydration free energy"
)
endpoint_order <- unname(endpoint_labels)
regime_labels <- c(
  "same assay, disjoint compounds" = "Same-assay source",
  "different ADME assay" = "Different-assay source"
)
method_labels <- c(
  Proposed = "Proposed",
  `F-AdaTrans` = "F-AdaTrans",
  `Target-only` = "Target-only",
  `Complete pooling` = "Complete pooling",
  TransLasso = "TransLasso"
)

tasks[, endpoint := endpoint_labels[drug]]
tasks[, regime_label := regime_labels[source_regime]]
task_summary <- merge(
  task_summary,
  tasks[, .(task, endpoint, regime_label)],
  by = "task",
  all.x = TRUE
)
task_summary[, endpoint := factor(endpoint, levels = endpoint_order)]
task_summary[, regime_label := factor(
  regime_label,
  levels = c("Same-assay source", "Different-assay source")
)]

cluster_bootstrap_interval <- function(values, B = 10000L) {
  n <- length(values)
  draws <- replicate(B, mean(sample(values, n, replace = TRUE)))
  c(
    estimate = mean(values),
    lower_95 = unname(quantile(draws, 0.025)),
    upper_95 = unname(quantile(draws, 0.975))
  )
}
set.seed(20261114L)
same_assay_contrast <- cluster_bootstrap_interval(
  task_summary[regime_label == "Same-assay source", Proposed - Target_only]
)
different_assay_contrast <- cluster_bootstrap_interval(
  task_summary[
    regime_label == "Different-assay source",
    Proposed - F_AdaTrans
  ]
)
endpoint_contrasts <- task_summary[, .(
  proposed_minus_target = mean(Proposed - Target_only),
  proposed_minus_transfer = mean(Proposed - F_AdaTrans),
  proposed_minus_translasso = mean(Proposed - TransLasso)
), by = endpoint]
overall_target_contrast <- cluster_bootstrap_interval(
  endpoint_contrasts$proposed_minus_target
)
overall_transfer_contrast <- cluster_bootstrap_interval(
  endpoint_contrasts$proposed_minus_transfer
)
overall_translasso_contrast <- cluster_bootstrap_interval(
  endpoint_contrasts$proposed_minus_translasso
)
regime_contrast_intervals <- data.table(
  comparison = c(
    "Same-assay: Proposed minus Target-only",
    "Different-assay: Proposed minus F-AdaTrans"
  ),
  estimate = c(same_assay_contrast["estimate"],
               different_assay_contrast["estimate"]),
  lower_95 = c(same_assay_contrast["lower_95"],
               different_assay_contrast["lower_95"]),
  upper_95 = c(same_assay_contrast["upper_95"],
               different_assay_contrast["upper_95"])
)
fwrite(
  regime_contrast_intervals,
  file.path(result_dir, "adme_nfit40_10rep_regime_contrast_intervals.csv")
)
overall_contrast_intervals <- data.table(
  comparison = c(
    "Proposed minus Target-only",
    "Proposed minus F-AdaTrans",
    "Proposed minus TransLasso"
  ),
  estimate = c(
    overall_target_contrast["estimate"],
    overall_transfer_contrast["estimate"],
    overall_translasso_contrast["estimate"]
  ),
  lower_95 = c(
    overall_target_contrast["lower_95"],
    overall_transfer_contrast["lower_95"],
    overall_translasso_contrast["lower_95"]
  ),
  upper_95 = c(
    overall_target_contrast["upper_95"],
    overall_transfer_contrast["upper_95"],
    overall_translasso_contrast["upper_95"]
  )
)
fwrite(
  overall_contrast_intervals,
  file.path(result_dir, "adme_nfit40_10rep_overall_contrast_intervals.csv")
)
task_summary[, selection_se := sqrt(
  selection_rate * (1 - selection_rate) / 10
)]
task_summary[, `:=`(
  selection_lower = pmax(0, selection_rate - selection_se),
  selection_upper = pmin(1, selection_rate + selection_se)
)]

base_theme <- theme_bw(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "bottom"
  )

p_selection <- ggplot(
  task_summary,
  aes(endpoint, selection_rate, colour = regime_label, group = regime_label)
) +
  geom_hline(yintercept = 0.5, linetype = 3, colour = "grey55") +
  geom_ribbon(
    aes(
      ymin = selection_lower,
      ymax = selection_upper,
      fill = regime_label
    ),
    colour = NA,
    alpha = 0.14
  ) +
  geom_line(linewidth = 0.45) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_colour_manual(values = c("#2166ac", "#b2182b"), name = NULL) +
  scale_fill_manual(
    values = c("#2166ac", "#b2182b"),
    name = NULL,
    guide = "none"
  ) +
  labs(x = NULL, y = "Transfer-selection rate") +
  base_theme
ggsave(
  file.path(figure_dir, "adme_transfer_selection.pdf"),
  p_selection,
  width = 7.4,
  height = 3.5
)

replication_wide <- dcast(
  main[method %chin% c("Proposed", "F-AdaTrans", "Target-only")],
  task + repetition ~ method,
  value.var = "prediction_mae"
)
replication_wide[, `:=`(
  difference_target = Proposed - `Target-only`,
  difference_transfer = Proposed - `F-AdaTrans`
)]
replication_wide <- merge(
  replication_wide,
  tasks[, .(task, endpoint, regime_label)],
  by = "task",
  all.x = TRUE
)
contrasts <- rbindlist(list(
  replication_wide[regime_label == "Same-assay source", .(
    difference = mean(difference_target),
    difference_sd = sd(difference_target)
  ), by = .(task, endpoint)][, comparison :=
    "Same-assay source: Proposed minus Target-only"],
  replication_wide[regime_label == "Different-assay source", .(
    difference = mean(difference_transfer),
    difference_sd = sd(difference_transfer)
  ), by = .(task, endpoint)][, comparison :=
    "Different-assay source: Proposed minus F-AdaTrans"]
))
contrasts[, endpoint := factor(endpoint, levels = endpoint_order)]
contrasts[, comparison := factor(
  comparison,
  levels = c(
    "Same-assay source: Proposed minus Target-only",
    "Different-assay source: Proposed minus F-AdaTrans"
  )
)]
p_contrast <- ggplot(
  contrasts,
  aes(endpoint, difference, group = 1)
) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey45") +
  geom_ribbon(
    aes(
      ymin = difference - difference_sd,
      ymax = difference + difference_sd
    ),
    fill = "#9ecae1",
    alpha = 0.35
  ) +
  geom_line(linewidth = 0.75, colour = "#2166ac") +
  facet_wrap(~comparison, ncol = 1) +
  labs(x = NULL, y = "MAE difference: Proposed minus comparator") +
  base_theme
ggsave(
  file.path(figure_dir, "adme_task_mae_contrasts.pdf"),
  p_contrast,
  width = 7.4,
  height = 5.7
)

interval_wide <- dcast(
  main[method %chin% c("Proposed", "F-AdaTrans", "Target-only")],
  task + repetition ~ method,
  value.var = "interval_length"
)
interval_wide[, `:=`(
  difference_target = Proposed - `Target-only`,
  difference_transfer = Proposed - `F-AdaTrans`
)]
interval_wide <- merge(
  interval_wide,
  tasks[, .(task, endpoint, regime_label)],
  by = "task",
  all.x = TRUE
)
interval_contrasts <- rbindlist(list(
  interval_wide[regime_label == "Same-assay source", .(
    difference = mean(difference_target),
    difference_sd = sd(difference_target)
  ), by = .(task, endpoint)][, comparison :=
    "Same-assay source: Proposed minus Target-only"],
  interval_wide[regime_label == "Different-assay source", .(
    difference = mean(difference_transfer),
    difference_sd = sd(difference_transfer)
  ), by = .(task, endpoint)][, comparison :=
    "Different-assay source: Proposed minus F-AdaTrans"]
))
interval_contrasts[, endpoint := factor(endpoint, levels = endpoint_order)]
interval_contrasts[, comparison := factor(
  comparison,
  levels = c(
    "Same-assay source: Proposed minus Target-only",
    "Different-assay source: Proposed minus F-AdaTrans"
  )
)]
p_interval_contrast <- ggplot(
  interval_contrasts,
  aes(endpoint, difference, group = 1)
) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey45") +
  geom_ribbon(
    aes(
      ymin = difference - difference_sd,
      ymax = difference + difference_sd
    ),
    fill = "#9ecae1",
    alpha = 0.35
  ) +
  geom_line(linewidth = 0.75, colour = "#2166ac") +
  facet_wrap(~comparison, ncol = 1) +
  labs(x = NULL, y = "Length difference: Proposed minus comparator") +
  base_theme
ggsave(
  file.path(figure_dir, "adme_task_interval_contrasts.pdf"),
  p_interval_contrast,
  width = 7.4,
  height = 5.7
)

results <- merge(
  results,
  tasks[, .(task, source_regime)],
  by = "task",
  all.x = TRUE
)
results[, regime_label := factor(
  regime_labels[source_regime],
  levels = c("Same-assay source", "Different-assay source")
)]
regime_within_task_sd <- results[
  , .(mae_sd = sd(prediction_mae)),
  by = .(source_regime, method, task)
]
regime_stability <- regime_within_task_sd[
  , .(mean_within_task_mae_sd = mean(mae_sd)),
  by = .(source_regime, method)
]
selection_stability <- task_summary[
  , .(
    selection_rate = mean(selection_rate),
    across_endpoint_sd = sd(selection_rate)
  ),
  by = regime_label
]
fwrite(
  regime_stability,
  file.path(result_dir, "adme_nfit40_10rep_regime_stability.csv")
)
fwrite(
  selection_stability,
  file.path(result_dir, "adme_nfit40_10rep_selection_stability.csv")
)
results[, method := factor(method, levels = names(method_labels))]
p_coverage <- ggplot(results, aes(method, coverage, fill = method)) +
  geom_hline(yintercept = 0.9, linetype = 2, colour = "grey35") +
  geom_boxplot(width = 0.68, outlier.size = 0.5) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  facet_wrap(~regime_label, nrow = 1) +
  labs(x = NULL, y = "Target-domain coverage", title = "(a) Coverage") +
  base_theme
p_length <- ggplot(results, aes(method, interval_length, fill = method)) +
  geom_boxplot(width = 0.68, outlier.size = 0.5) +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  facet_wrap(~regime_label, nrow = 1) +
  labs(x = NULL, y = "Prediction-interval length", title = "(b) Length") +
  base_theme
pdf(
  file.path(figure_dir, "adme_coverage_length_distributions.pdf"),
  width = 11.2,
  height = 4.2
)
grid.arrange(p_coverage, p_length, nrow = 1)
dev.off()
