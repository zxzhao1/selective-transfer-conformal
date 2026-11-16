#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

root <- normalizePath(getwd(), mustWork = TRUE)
result_dir <- file.path(root, "results")
main_file <- file.path(
  result_dir, "real_results_adme_nfit40_10rep_p300.csv"
)
translasso_file <- file.path(
  result_dir, "real_results_translasso_adme_nfit40_10rep_p300.csv"
)
task_file <- file.path(root, "data", "processed", "adme_tasks.csv")
stopifnot(file.exists(main_file), file.exists(translasso_file),
          file.exists(task_file))

main <- fread(main_file)
translasso <- fread(translasso_file)
tasks <- fread(task_file)
results <- rbind(main, translasso, fill = TRUE)

stopifnot(
  nrow(main) == 18L * 10L * 4L,
  nrow(translasso) == 18L * 10L,
  nrow(results) == 18L * 10L * 5L,
  uniqueN(main[, .(task, repetition, method)]) == nrow(main),
  uniqueN(translasso[, .(task, repetition, method)]) == nrow(translasso),
  setequal(unique(main$method), c(
    "Proposed", "F-AdaTrans", "Target-only", "Complete pooling"
  )),
  identical(unique(translasso$method), "TransLasso"),
  all(results$coverage >= 0 & results$coverage <= 1),
  all(is.finite(results$interval_length) & results$interval_length >= 0),
  all(is.finite(results$prediction_mae)),
  all(is.finite(results$prediction_mse))
)

results <- merge(
  results,
  tasks[, .(task, drug, source_regime)],
  by = "task",
  all.x = TRUE,
  sort = FALSE
)
stopifnot(!anyNA(results$source_regime))

metric_names <- c(
  "coverage", "interval_length", "prediction_mse", "prediction_mae"
)
method_summary <- results[, c(
  list(task_replications = .N),
  lapply(.SD, mean),
  lapply(.SD, sd)
), by = method, .SDcols = metric_names]
setnames(
  method_summary,
  3L:ncol(method_summary),
  c(paste0("mean_", metric_names), paste0("sd_", metric_names))
)
method_order <- c(
  "Proposed", "F-AdaTrans", "Target-only", "Complete pooling", "TransLasso"
)
method_summary[, method_rank := match(method, method_order)]
setorder(method_summary, method_rank)
method_summary[, method_rank := NULL]
fwrite(
  method_summary,
  file.path(result_dir, "adme_nfit40_10rep_method_summary.csv")
)

selection <- results[method == "Proposed", .(
  selection_rate = mean(selected_transfer),
  selection_rate_sd = sd(selected_transfer)
), by = .(task, source_regime)]
task_metrics <- results[, lapply(.SD, mean),
  by = .(task, source_regime, method), .SDcols = metric_names]

task_mae <- dcast(
  task_metrics, task + source_regime ~ method,
  value.var = "prediction_mae"
)
setnames(
  task_mae,
  c("F-AdaTrans", "Target-only", "Complete pooling"),
  c("F_AdaTrans", "Target_only", "Complete_pooling")
)
task_summary <- merge(task_mae, selection, by = c("task", "source_regime"))
setcolorder(task_summary, c(
  "task", "source_regime", "selection_rate", "selection_rate_sd",
  "Proposed", "F_AdaTrans", "Target_only", "Complete_pooling", "TransLasso"
))
setorder(task_summary, task)
fwrite(
  task_summary,
  file.path(result_dir, "adme_nfit40_10rep_task_summary.csv")
)

regime_metrics <- results[, lapply(.SD, mean),
  by = .(source_regime, method), .SDcols = metric_names]
regime_selection <- selection[, .(
  selection_rate = mean(selection_rate),
  across_task_sd = sd(selection_rate)
), by = source_regime]
regime_summary <- merge(
  regime_metrics,
  regime_selection,
  by = "source_regime",
  all.x = TRUE
)
regime_summary[, method_rank := match(method, method_order)]
setorder(regime_summary, source_regime, method_rank)
regime_summary[, method_rank := NULL]
fwrite(
  regime_summary,
  file.path(result_dir, "adme_nfit40_10rep_regime_summary.csv")
)

paired <- dcast(
  results[method %chin% c("Proposed", "F-AdaTrans", "Target-only")],
  task + drug + source_regime + repetition ~ method,
  value.var = c("coverage", "interval_length", "prediction_mae")
)
paired[, `:=`(
  interval_oracle = pmin(`interval_length_F-AdaTrans`,
                         `interval_length_Target-only`),
  interval_regret = interval_length_Proposed -
    pmin(`interval_length_F-AdaTrans`, `interval_length_Target-only`),
  length_minus_transfer = interval_length_Proposed -
    `interval_length_F-AdaTrans`,
  length_minus_target = interval_length_Proposed - `interval_length_Target-only`,
  coverage_minus_transfer = coverage_Proposed - `coverage_F-AdaTrans`,
  coverage_minus_target = coverage_Proposed - `coverage_Target-only`,
  mae_minus_transfer = prediction_mae_Proposed - `prediction_mae_F-AdaTrans`,
  mae_minus_target = prediction_mae_Proposed - `prediction_mae_Target-only`
)]
fwrite(
  paired,
  file.path(result_dir, "adme_nfit40_10rep_paired_metrics.csv")
)

regime_interval_summary <- paired[, .(
  repetitions = .N,
  proposed_coverage_mean = mean(coverage_Proposed),
  proposed_coverage_sd = sd(coverage_Proposed),
  proposed_length_mean = mean(interval_length_Proposed),
  proposed_length_sd = sd(interval_length_Proposed),
  interval_regret_mean = mean(interval_regret),
  interval_regret_sd = sd(interval_regret),
  length_minus_transfer_mean = mean(length_minus_transfer),
  length_minus_target_mean = mean(length_minus_target)
), by = source_regime]
fwrite(
  regime_interval_summary,
  file.path(result_dir, "adme_nfit40_10rep_interval_summary.csv")
)

bootstrap_mean <- function(x, B = 10000L) {
  n <- length(x)
  draws <- replicate(B, mean(sample(x, n, replace = TRUE)))
  data.table(
    estimate = mean(x),
    lower_95 = unname(quantile(draws, 0.025)),
    upper_95 = unname(quantile(draws, 0.975))
  )
}

interval_regime_contrasts <- rbindlist(list(
  paired[source_regime == "same assay, disjoint compounds", .(
    value = mean(length_minus_target)
  ), by = drug][, comparison :=
    "Same-assay: Proposed minus Target-only"],
  paired[source_regime == "different ADME assay", .(
    value = mean(length_minus_transfer)
  ), by = drug][, comparison :=
    "Different-assay: Proposed minus F-AdaTrans"]
))
set.seed(20261115L)
interval_regime_intervals <- interval_regime_contrasts[, {
  interval <- bootstrap_mean(value)
  interval
}, by = comparison]
fwrite(
  interval_regime_intervals,
  file.path(result_dir, "adme_nfit40_10rep_interval_contrast_intervals.csv")
)

endpoint_contrasts <- paired[, .(
  proposed_minus_target = mean(prediction_mae_Proposed -
                                 `prediction_mae_Target-only`),
  proposed_minus_transfer = mean(prediction_mae_Proposed -
                                   `prediction_mae_F-AdaTrans`),
  interval_regret = mean(interval_regret)
), by = drug]
set.seed(20261115L)
contrast_list <- list(
  proposed_minus_target = endpoint_contrasts$proposed_minus_target,
  proposed_minus_transfer = endpoint_contrasts$proposed_minus_transfer,
  interval_regret = endpoint_contrasts$interval_regret
)
task_clustered_intervals <- rbindlist(lapply(
  names(contrast_list),
  function(name) bootstrap_mean(contrast_list[[name]])[, contrast := name]
), use.names = TRUE)
setcolorder(task_clustered_intervals, c(
  "contrast", "estimate", "lower_95", "upper_95"
))
fwrite(
  task_clustered_intervals,
  file.path(result_dir, "adme_nfit40_10rep_task_clustered_intervals.csv")
)

cat("ADME summaries validated and written.\n")
