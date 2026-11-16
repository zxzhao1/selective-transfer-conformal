#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
repetitions <- if (length(args) >= 1L) as.integer(args[[1L]]) else 10L
p_keep <- if (length(args) >= 2L) as.integer(args[[2L]]) else 128L
target_fit_cap <- if (length(args) >= 3L) as.integer(args[[3L]]) else 40L
author_file <- if (length(args) >= 4L) args[[4L]] else {
  file.path("external", "TransLasso", "TransLasso-functions.R")
}
design <- if (length(args) >= 5L) args[[5L]] else "gas_sensor"
stopifnot(design %in% c(
  "gas_sensor", "acs", "parkinsons", "adme", "mice_protein", "kiba"
))

root <- normalizePath(getwd(), mustWork = TRUE)
author_file <- normalizePath(author_file, mustWork = TRUE)

pipeline <- new.env(parent = globalenv())
pipeline$commandArgs <- function(trailingOnly = FALSE) {
  if (trailingOnly) {
    c(as.character(repetitions), as.character(p_keep), design,
      as.character(target_fit_cap))
  } else {
    commandArgs(trailingOnly = FALSE)
  }
}
old_definitions_only <- Sys.getenv("REAL_DATA_DEFINITIONS_ONLY", unset = NA)
Sys.setenv(REAL_DATA_DEFINITIONS_ONLY = "true")
on.exit({
  if (is.na(old_definitions_only)) {
    Sys.unsetenv("REAL_DATA_DEFINITIONS_ONLY")
  } else {
    Sys.setenv(REAL_DATA_DEFINITIONS_ONLY = old_definitions_only)
  }
}, add = TRUE)
sys.source(file.path(root, "02_run_pilot.R"), envir = pipeline)
source(file.path(dirname(root), "simulation", "06_translasso.R"))

tasks <- pipeline$tasks
output_file <- file.path(
  root, "results",
  sprintf(
    "real_results_translasso_%s_nfit%d_%drep_p%d.csv",
    design, target_fit_cap, repetitions, p_keep
  )
)
results <- if (file.exists(output_file)) fread(output_file) else data.table()

for (task_index in seq_len(nrow(tasks))) {
  spec <- tasks[task_index]
  for (rep in seq_len(repetitions)) {
    if (nrow(results) > 0L &&
        results[task == spec$task & repetition == rep, .N] == 1L) {
      next
    }
    seed <- 20261114L + 1000L * task_index + rep
    cat("Task", spec$task, "replication", rep, "of", repetitions, "\n")
    data <- pipeline$prepare_replication(spec, seed, p_keep)
    cfg <- pipeline$make_real_config(data, seed)
    fitted <- fit_and_calibrate_translasso(
      data, cfg, author_file, seed + 500L
    )
    row <- pipeline$evaluate(
      fitted, data$target_test, data$metric_scale, spec$task, rep
    )
    if (nrow(results) > 0L) {
      results <- results[!(task == spec$task & repetition == rep)]
    }
    results <- rbind(results, row, fill = TRUE)
    setorder(results, task, repetition)
    fwrite(results, output_file)
  }
}

stopifnot(
  nrow(results) == nrow(tasks) * repetitions,
  results[, uniqueN(task)] == nrow(tasks),
  results[, uniqueN(repetition)] == repetitions,
  results[, uniqueN(method)] == 1L,
  results[, unique(method)] == "TransLasso"
)
print(results)
