#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

options(error = function() {
  traceback(30L)
  quit(save = "no", status = 1L)
})

args <- commandArgs(trailingOnly = TRUE)
repetitions <- if (length(args) >= 1L) as.integer(args[1L]) else 1L
p_keep <- if (length(args) >= 2L) as.integer(args[2L]) else 100L
design_override <- Sys.getenv("REAL_DATA_DESIGN", unset = "")
fit_cap_override <- Sys.getenv("REAL_DATA_TARGET_FIT_CAP", unset = "")
design <- if (nzchar(design_override)) {
  design_override
} else if (length(args) >= 3L) {
  args[3L]
} else {
  "cancer"
}
target_fit_cap <- if (nzchar(fit_cap_override)) {
  as.integer(fit_cap_override)
} else if (length(args) >= 4L) {
  as.integer(args[4L])
} else {
  NA_integer_
}
stopifnot(design %chin% c(
  "cancer", "system", "system_panel", "measurement"
  , "heterogeneity", "cross_study", "beataml", "beataml_ctrp",
  "gas_sensor", "acs", "parkinsons", "adme", "mice_protein", "kiba"
))

root <- normalizePath(getwd(), mustWork = TRUE)
simulation_dir <- file.path(dirname(root), "simulation")
for (file in c(
  "00_config.R", "01_utils.R", "03_estimators.R", "04_proposed.R"
)) {
  source(file.path(simulation_dir, file))
}
require_simulation_packages()

expression_file <- if (design == "cross_study") {
  "cross_study_expression.rds"
} else if (design == "beataml") {
  "beataml_expression.rds"
} else if (design == "beataml_ctrp") {
  "beataml_ctrp_expression.rds"
} else if (design == "gas_sensor") {
  "gas_sensor_expression.rds"
} else if (design == "acs") {
  "acs_expression.rds"
} else if (design == "parkinsons") {
  "parkinsons_expression.rds"
} else if (design == "adme") {
  "adme_expression.rds"
} else if (design == "mice_protein") {
  "mice_protein_expression.rds"
} else if (design == "kiba") {
  "kiba_expression.rds"
} else {
  "gdsc_expression.rds"
}
expression <- readRDS(file.path(root, "data", "processed", expression_file))
if (design == "kiba") {
  response <- fread(file.path(
    root, "data", "processed", "kiba_response.csv"
  ))
  tasks <- fread(file.path(
    root, "data", "processed", "kiba_tasks.csv"
  ))
  domain_variable <- "domain"
} else if (design == "mice_protein") {
  response <- fread(file.path(
    root, "data", "processed", "mice_protein_response.csv"
  ))
  tasks <- fread(file.path(
    root, "data", "processed", "mice_protein_tasks.csv"
  ))
  domain_variable <- "domain"
} else if (design == "adme") {
  response <- fread(file.path(
    root, "data", "processed", "adme_response.csv"
  ))
  tasks <- fread(file.path(
    root, "data", "processed", "adme_tasks.csv"
  ))
  domain_variable <- "domain"
} else if (design == "parkinsons") {
  response <- fread(file.path(
    root, "data", "processed", "parkinsons_response.csv"
  ))
  tasks <- fread(file.path(
    root, "data", "processed", "parkinsons_tasks.csv"
  ))
  domain_variable <- "domain"
} else if (design == "acs") {
  response <- fread(file.path(
    root, "data", "processed", "acs_response.csv"
  ))
  tasks <- fread(file.path(
    root, "data", "processed", "acs_tasks.csv"
  ))
  domain_variable <- "domain"
} else if (design == "gas_sensor") {
  response <- fread(file.path(
    root, "data", "processed", "gas_sensor_response.csv"
  ))
  tasks <- fread(file.path(
    root, "data", "processed", "gas_sensor_tasks.csv"
  ))
  domain_variable <- "domain"
} else if (design == "beataml") {
  response <- fread(file.path(
    root, "data", "processed", "beataml_response.csv"
  ))
  tasks <- fread(file.path(
    root, "data", "processed", "beataml_tasks.csv"
  ))
  domain_variable <- "domain"
} else if (design == "beataml_ctrp") {
  response <- fread(file.path(
    root, "data", "processed", "beataml_ctrp_response.csv"
  ))
  tasks <- fread(file.path(
    root, "data", "processed", "beataml_ctrp_tasks.csv"
  ))
  domain_variable <- "domain"
} else if (design == "cross_study") {
  response <- fread(file.path(
    root, "data", "processed", "cross_study_response.csv"
  ))
  tasks <- fread(file.path(
    root, "data", "processed", "cross_study_tasks.csv"
  ))
  domain_variable <- "domain"
} else if (design == "heterogeneity") {
  response <- fread(file.path(
    root, "data", "processed", "gdsc_heterogeneity_response.csv"
  ))
  tasks <- fread(file.path(
    root, "data", "processed", "gdsc_heterogeneity_tasks.csv"
  ))
  domain_variable <- "domain"
} else if (design == "measurement") {
  response <- fread(file.path(
    root, "data", "processed", "gdsc_measurement_response.csv"
  ))
  tasks <- fread(file.path(
    root, "data", "processed", "gdsc_measurement_tasks.csv"
  ))
  domain_variable <- "analysis_domain"
} else if (design %chin% c("system", "system_panel")) {
  response <- fread(file.path(
    root, "data", "processed", "gdsc_system_response.csv"
  ))
  task_file <- if (design == "system_panel") {
    "gdsc_system_panel_tasks.csv"
  } else {
    "gdsc_system_tasks.csv"
  }
  tasks <- fread(file.path(root, "data", "processed", task_file))
  domain_variable <- "domain"
} else {
  response <- fread(file.path(root, "data", "processed", "gdsc_response.csv"))
  tasks <- fread(file.path(root, "data", "processed", "gdsc_tasks.csv"))
  domain_variable <- "cancer_type"
}

make_split <- function(
    n, seed, groups = NULL, fit_cap = NA_integer_, strata = NULL
) {
  set.seed(seed)
  if (!is.null(groups)) {
    stopifnot(length(groups) == n, !anyNA(groups))
    group_rows <- split(seq_len(n), as.character(groups))
    group_order <- sample(names(group_rows))
    group_sizes <- lengths(group_rows)
    desired_fit <- if (is.finite(fit_cap)) {
      min(as.integer(fit_cap), floor(0.60 * n))
    } else {
      floor(0.60 * n)
    }
    desired_fit <- max(20L, desired_fit)
    choose_count <- function(order, desired, max_count) {
      candidates <- seq_len(max_count)
      cumulative <- cumsum(group_sizes[order])
      candidates[which.min(abs(cumulative[candidates] - desired))]
    }
    fit_count <- choose_count(
      group_order, desired_fit, max(1L, length(group_order) - 2L)
    )
    fit_groups <- group_order[seq_len(fit_count)]
    remaining <- setdiff(group_order, fit_groups)
    desired_cal <- max(9L, floor(0.20 * n))
    cal_count <- choose_count(
      remaining, desired_cal, max(1L, length(remaining) - 1L)
    )
    cal_groups <- remaining[seq_len(cal_count)]
    test_groups <- setdiff(remaining, cal_groups)
    return(list(
      fit = sort(unlist(group_rows[fit_groups], use.names = FALSE)),
      cal = sort(unlist(group_rows[cal_groups], use.names = FALSE)),
      test = sort(unlist(group_rows[test_groups], use.names = FALSE))
    ))
  }
  n_fit <- if (is.finite(fit_cap)) {
    min(as.integer(fit_cap), floor(0.60 * n))
  } else {
    floor(0.60 * n)
  }
  n_fit <- max(20L, n_fit)
  n_cal <- max(9L, floor(0.20 * n))
  if (n_fit + n_cal >= n) {
    n_fit <- floor(n / 3)
    n_cal <- floor(n / 3)
  }
  if (!is.null(strata)) {
    stopifnot(length(strata) == n, !anyNA(strata))
    partition <- make_stratified_partition(
      strata, c(n_fit, n_cal, n - n_fit - n_cal)
    )
    shuffle <- function(rows) rows[sample.int(length(rows))]
    index <- c(
      shuffle(which(partition == 1L)),
      shuffle(which(partition == 2L)),
      shuffle(which(partition == 3L))
    )
  } else {
    index <- sample.int(n)
  }
  list(
    fit = index[seq_len(n_fit)],
    cal = index[n_fit + seq_len(n_cal)],
    test = index[(n_fit + n_cal + 1L):n]
  )
}

prepare_replication <- function(spec, seed, p_keep) {
  source_types <- strsplit(spec$sources, ",", fixed = TRUE)[[1L]]
  domains <- c(spec$target, source_types)
  dat <- response[
    drug == spec$drug & get(domain_variable) %chin% domains &
      expression_id %chin% rownames(expression)
  ]
  setorderv(dat, c(domain_variable, "expression_id"))
  dat <- unique(dat, by = c("drug", domain_variable, "expression_id"))

  target_rows <- dat[get(domain_variable) == spec$target]
  split_groups <- if (design %chin% c(
    "beataml", "beataml_ctrp", "mice_protein"
  )) {
    target_rows$group_id
  } else if (design == "adme") {
    stopifnot("compound_id" %chin% names(target_rows))
    target_rows$compound_id
  } else {
    NULL
  }
  split <- make_split(
    nrow(target_rows), seed, groups = split_groups,
    fit_cap = target_fit_cap,
    strata = if (design == "gas_sensor") target_rows$response_value else NULL
  )
  if (design == "adme") {
    fit_compounds <- target_rows$compound_id[split$fit]
    cal_compounds <- target_rows$compound_id[split$cal]
    test_compounds <- target_rows$compound_id[split$test]
    stopifnot(
      !anyDuplicated(fit_compounds),
      !anyDuplicated(cal_compounds),
      !anyDuplicated(test_compounds),
      length(intersect(fit_compounds, cal_compounds)) == 0L,
      length(intersect(fit_compounds, test_compounds)) == 0L,
      length(intersect(cal_compounds, test_compounds)) == 0L
    )
  }
  source_rows <- lapply(
    source_types, function(type) dat[get(domain_variable) == type]
  )
  if (design == "cross_study") {
    source_rows <- lapply(source_rows, function(rows) {
      rows[!expression_id %chin% target_rows$expression_id]
    })
  } else if (design == "adme") {
    stopifnot("compound_id" %chin% names(target_rows))
    target_compounds <- unique(target_rows$compound_id)
    source_rows <- lapply(source_rows, function(rows) {
      stopifnot("compound_id" %chin% names(rows))
      rows[!compound_id %chin% target_compounds]
    })
    stopifnot(all(vapply(
      source_rows,
      function(rows) !any(rows$compound_id %chin% target_compounds),
      logical(1)
    )))
  }

  response_column <- if ("response_value" %chin% names(dat)) {
    "response_value"
  } else {
    "ln_ic50"
  }

  domain_standardization <- design %chin% c(
    "beataml", "beataml_ctrp", "gas_sensor", "acs", "parkinsons", "adme",
    "mice_protein", "kiba"
  )
  training_ids <- if (domain_standardization) {
    target_rows$expression_id[split$fit]
  } else {
    c(
      target_rows$expression_id[split$fit],
      unlist(lapply(source_rows, `[[`, "expression_id"), use.names = FALSE)
    )
  }
  training_x <- expression[training_ids, , drop = FALSE]
  variance <- apply(training_x, 2L, stats::var, na.rm = TRUE)
  eligible <- is.finite(variance) & variance > 1e-8
  selected <- names(sort(variance[eligible], decreasing = TRUE))[
    seq_len(min(p_keep, sum(eligible)))
  ]

  y_center <- mean(target_rows[[response_column]][split$fit])
  y_scale <- stats::sd(target_rows[[response_column]][split$fit])
  if (!is.finite(y_scale) || y_scale <= 0) stop("Invalid target response scale.")
  y_for <- function(y) (y - y_center) / y_scale

  if (domain_standardization) {
    target_fit_x <- expression[
      target_rows$expression_id[split$fit], selected, drop = FALSE
    ]
    target_mean <- colMeans(target_fit_x)
    target_sd <- apply(target_fit_x, 2L, stats::sd)
    target_sd[!is.finite(target_sd) | target_sd <= 1e-8] <- 1
    standardize_x <- function(x, center, scale) {
      x <- sweep(x, 2L, center, "-")
      x <- sweep(x, 2L, scale, "/")
      x[!is.finite(x)] <- 0
      x
    }
    target_fit <- list(
      x = standardize_x(target_fit_x, target_mean, target_sd),
      y = y_for(target_rows[[response_column]][split$fit])
    )
    target_cal <- list(
      x = standardize_x(expression[
        target_rows$expression_id[split$cal], selected, drop = FALSE
      ], target_mean, target_sd),
      y = y_for(target_rows[[response_column]][split$cal])
    )
    target_test <- list(
      x = standardize_x(expression[
        target_rows$expression_id[split$test], selected, drop = FALSE
      ], target_mean, target_sd),
      y = y_for(target_rows[[response_column]][split$test])
    )
    sources <- lapply(source_rows, function(rows) {
      source_x <- expression[rows$expression_id, selected, drop = FALSE]
      if (design == "adme") {
        source_mean <- target_mean
        source_sd <- target_sd
      } else {
        source_mean <- colMeans(source_x)
        source_sd <- apply(source_x, 2L, stats::sd)
        source_sd[!is.finite(source_sd) | source_sd <= 1e-8] <- 1
      }
      same_adme_assay <- design == "adme" &&
        identical(as.character(spec$source_regime),
                  "same assay, disjoint compounds")
      source_y_mean <- if (same_adme_assay) {
        y_center
      } else {
        mean(rows[[response_column]])
      }
      source_y_sd <- if (same_adme_assay) {
        y_scale
      } else {
        stats::sd(rows[[response_column]])
      }
      if (!is.finite(source_y_sd) || source_y_sd <= 0) {
        stop("Invalid source response scale.")
      }
      list(
        x = standardize_x(source_x, source_mean, source_sd),
        y = (rows[[response_column]] - source_y_mean) / source_y_sd
      )
    })
  } else {
    x_for <- function(ids) expression[ids, selected, drop = FALSE]
    target_fit <- list(
      x = x_for(target_rows$expression_id[split$fit]),
      y = y_for(target_rows[[response_column]][split$fit])
    )
    target_cal <- list(
      x = x_for(target_rows$expression_id[split$cal]),
      y = y_for(target_rows[[response_column]][split$cal])
    )
    target_test <- list(
      x = x_for(target_rows$expression_id[split$test]),
      y = y_for(target_rows[[response_column]][split$test])
    )
    sources <- lapply(source_rows, function(rows) list(
      x = x_for(rows$expression_id),
      y = y_for(rows[[response_column]])
    ))
  }

  list(
    target_fit = target_fit,
    target_cal = target_cal,
    target_test = target_test,
    sources = sources,
    y_center = y_center,
    y_scale = y_scale,
    metric_scale = if (domain_standardization) 1 else y_scale,
    selected_genes = selected,
    counts = c(
      target_fit = nrow(target_fit$x),
      target_cal = nrow(target_cal$x),
      target_test = nrow(target_test$x),
      vapply(sources, function(z) nrow(z$x), integer(1))
    )
  )
}

make_real_config <- function(data, seed) {
  cfg <- make_config("smoke")
  cfg$p <- ncol(data$target_fit$x)
  cfg$K <- length(data$sources)
  cfg$n_fit <- nrow(data$target_fit$x)
  cfg$n_cal <- nrow(data$target_cal$x)
  cfg$n_test <- nrow(data$target_test$x)
  cfg$n_source <- vapply(data$sources, function(z) nrow(z$x), integer(1))
  cfg$outer_folds <- 3L
  cfg$inner_folds <- 3L
  cfg$lambda_multipliers <- c(0.7, 1.4)
  cfg$linked_transfer_penalties <- TRUE
  cfg$initial_tuning <- "joint"
  cfg$scad_a <- max(3.7, 3.7 * cfg$K / 2)
  cfg$seed <- seed
  validate_config(cfg)
}

evaluate <- function(calibrated, test, y_scale, task, repetition) {
  prediction <- predict_linear_model(calibrated$model, test$x)
  lower <- prediction - calibrated$q
  upper <- prediction + calibrated$q
  data.table(
    task = task,
    repetition = repetition,
    method = calibrated$model$method,
    coverage = mean(test$y >= lower & test$y <= upper),
    interval_length = 2 * calibrated$q * y_scale,
    prediction_mse = mean((test$y - prediction)^2) * y_scale^2,
    prediction_mae = mean(abs(test$y - prediction)) * y_scale,
    selected_transfer = if (is.null(calibrated$selected)) NA_integer_ else
      as.integer(calibrated$selected == "transfer"),
    selection_risk_transfer = (calibrated$r_transfer %||% NA_real_) * y_scale,
    selection_risk_target = (calibrated$r_target %||% NA_real_) * y_scale,
    selection_risk_difference = (
      (calibrated$r_transfer %||% NA_real_) -
        (calibrated$r_target %||% NA_real_)
    ) * y_scale,
    lambda_initial = calibrated$model$lambda_initial %||% NA_real_,
    lambda0 = calibrated$model$lambda0 %||% NA_real_,
    lambda1 = calibrated$model$lambda1 %||% NA_real_,
    lambda_target = calibrated$model$lambda %||% NA_real_
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

run_main <- Sys.getenv("REAL_DATA_DEFINITIONS_ONLY") != "true"
if (run_main) {
dir.create(file.path(root, "results"), showWarnings = FALSE)
design_tag <- if (is.finite(target_fit_cap)) {
  sprintf("%s_nfit%d", design, target_fit_cap)
} else {
  design
}
result_checkpoint <- file.path(
  root, "results",
  sprintf("real_results_%s_%drep_p%d.csv", design_tag, repetitions, p_keep)
)
count_checkpoint <- file.path(
  root, "results",
  sprintf("real_sample_sizes_%s_%drep_p%d.csv", design_tag, repetitions, p_keep)
)
result_table <- if (file.exists(result_checkpoint)) {
  fread(result_checkpoint)
} else {
  data.table()
}
count_table <- if (file.exists(count_checkpoint)) {
  fread(count_checkpoint)
} else {
  data.table()
}

for (task_index in seq_len(nrow(tasks))) {
  spec <- tasks[task_index]
  for (rep in seq_len(repetitions)) {
    already_complete <- nrow(result_table) > 0L &&
      result_table[task == spec$task & repetition == rep, .N] == 4L &&
      nrow(count_table) > 0L &&
      count_table[task == spec$task & repetition == rep, .N] == 1L
    if (already_complete) {
      cat("Skipping completed task", spec$task, "replication", rep, "\n")
      next
    }
    seed <- 20261114L + 1000L * task_index + rep
    cat("Task", spec$task, "replication", rep, "of", repetitions, "\n")
    data <- prepare_replication(spec, seed, p_keep)
    cfg <- make_real_config(data, seed)
    source_initials <- fit_source_initials(data$sources, cfg, seed + 10L)
    proposed_seed <- seed + 20L
    shared_candidate_seed <- proposed_seed + 300L

    fitted <- list(
      Proposed = fit_proposed(
        data$target_fit, data$target_cal, data$sources, cfg,
        source_initials, proposed_seed
      ),
      `F-AdaTrans` = fit_and_calibrate(
        "F-AdaTrans", data, cfg, source_initials, shared_candidate_seed
      ),
      `Target-only` = fit_and_calibrate(
        "Target-only", data, cfg, source_initials, shared_candidate_seed
      ),
      `Complete pooling` = fit_and_calibrate(
        "Complete pooling", data, cfg, source_initials, seed + 50L
      )
    )

    selected_candidate <- if (fitted$Proposed$selected == "transfer") {
      fitted$`F-AdaTrans`
    } else {
      fitted$`Target-only`
    }
    stopifnot(isTRUE(all.equal(
      predict_linear_model(fitted$Proposed$model, data$target_test$x),
      predict_linear_model(selected_candidate$model, data$target_test$x),
      tolerance = 1e-10
    )))
    stopifnot(isTRUE(all.equal(
      fitted$Proposed$q, selected_candidate$q, tolerance = 1e-10
    )))

    replication_results <- rbindlist(lapply(names(fitted), function(method) {
      evaluate(
        fitted[[method]], data$target_test, data$metric_scale, spec$task, rep
      )
    }), fill = TRUE)
    replication_counts <- data.table(
      task = spec$task,
      repetition = rep,
      p = ncol(data$target_fit$x),
      n_fit = nrow(data$target_fit$x),
      n_cal = nrow(data$target_cal$x),
      n_test = nrow(data$target_test$x),
      source_sizes = paste(
        vapply(data$sources, function(z) nrow(z$x), integer(1)),
        collapse = ","
      )
    )
    if (nrow(result_table) > 0L) {
      result_table <- result_table[!(task == spec$task & repetition == rep)]
    }
    if (nrow(count_table) > 0L) {
      count_table <- count_table[!(task == spec$task & repetition == rep)]
    }
    result_table <- rbind(result_table, replication_results, fill = TRUE)
    count_table <- rbind(count_table, replication_counts, fill = TRUE)
    setorder(result_table, task, repetition, method)
    setorder(count_table, task, repetition)
    fwrite(result_table, result_checkpoint)
    fwrite(count_table, count_checkpoint)
  }
}

fwrite(result_table, file.path(root, "results", "pilot_results.csv"))
fwrite(count_table, file.path(root, "results", "pilot_sample_sizes.csv"))
print(result_table)
}

