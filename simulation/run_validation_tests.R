script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1L])))

for (file in c(
  "00_config.R", "01_utils.R", "02_dgp.R", "03_estimators.R",
  "04_proposed.R", "05_metrics.R", "06_translasso.R", "07_scenarios.R"
)) {
  source(file.path(script_dir, file), local = FALSE)
}

require_simulation_packages()
cfg <- make_config("smoke")
data <- generate_simulation_data(cfg, seed = cfg$seed + 700L)
source_initials <- fit_source_initials(data$sources, cfg, cfg$seed + 710L)

fit_original <- fit_proposed(
  data$target_fit, data$target_cal, data$sources, cfg,
  source_initials, cfg$seed + 720L
)

matching_method <- if (fit_original$selected == "transfer") {
  "F-AdaTrans"
} else {
  "Target-only"
}
matching_candidate <- fit_and_calibrate(
  matching_method,
  data,
  cfg,
  source_initials,
  cfg$seed + 1020L
)

altered_calibration <- data$target_cal
altered_calibration$y <- altered_calibration$y + 100
fit_altered <- fit_proposed(
  data$target_fit, altered_calibration, data$sources, cfg,
  source_initials, cfg$seed + 720L
)

stopifnot(
  identical(fit_original$selected, fit_altered$selected),
  isTRUE(all.equal(fit_original$model$beta, fit_altered$model$beta)),
  isTRUE(all.equal(fit_original$model$intercept, fit_altered$model$intercept)),
  isTRUE(all.equal(fit_original$r_transfer, fit_altered$r_transfer)),
  isTRUE(all.equal(fit_original$r_target, fit_altered$r_target)),
  !isTRUE(all.equal(fit_original$q, fit_altered$q)),
  isTRUE(all.equal(
    fit_original$model$beta,
    matching_candidate$model$beta
  )),
  isTRUE(all.equal(
    fit_original$model$intercept,
    matching_candidate$model$intercept
  )),
  isTRUE(all.equal(fit_original$q, matching_candidate$q))
)

diagnostic_row <- evaluate_calibrated_method(
  fit_original, data$target_test, data$truth$beta0
)
stopifnot(
  isTRUE(all.equal(
    diagnostic_row$selection_radius_transfer,
    fit_original$r_transfer
  )),
  isTRUE(all.equal(
    diagnostic_row$selection_radius_target,
    fit_original$r_target
  )),
  isTRUE(all.equal(
    diagnostic_row$selection_radius_margin,
    fit_original$r_target - fit_original$r_transfer
  ))
)

scores <- c(1, 2, 3, 4)
stopifnot(
  identical(finite_sample_quantile(scores, 0.5), 3),
  is.infinite(finite_sample_quantile(scores, 0.01))
)

zero_signal <- generate_simulation_data(
  cfg, seed = cfg$seed + 730L, h_delta = 0, s_delta = cfg$s_delta
)
stopifnot(all(vapply(
  zero_signal$truth$deltas,
  function(delta) all(delta == 0),
  logical(1)
)))

holdouts <- make_translasso_holdouts(cfg$n_fit, cfg$seed + 740L)
grid_check <- make_lambda_grids(
  cfg,
  n_target = cfg$n_fit,
  total_n = cfg$n_fit + cfg$K * cfg$n_source,
  n_source = cfg$n_source
)
joint_grid_check <- expand.grid(
  lambda_initial = grid_check$initial,
  lambda0 = grid_check$lambda0,
  lambda1 = grid_check$lambda1,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
expected_lambda1 <- cfg$lambda_multipliers *
  (cfg$n_source / (cfg$n_fit + cfg$K * cfg$n_source)) *
  sqrt(log(cfg$p) / cfg$n_fit)
full_cfg <- make_config("full")
simulation4 <- make_scenarios(4L, full_cfg)
stopifnot(
  length(holdouts$first) == floor(cfg$n_fit / 3L),
  length(holdouts$second) == floor(cfg$n_fit / 3L),
  length(intersect(holdouts$first, holdouts$second)) == 0L,
  all(c(holdouts$first, holdouts$second) >= 1L),
  all(c(holdouts$first, holdouts$second) <= cfg$n_fit),
  isTRUE(all.equal(grid_check$lambda1, expected_lambda1)),
  !isTRUE(cfg$linked_transfer_penalties),
  nrow(joint_grid_check) == length(cfg$lambda_multipliers)^3,
  "Ora-F-Ada" %in% methods_for_scenario(4L, simulation4[[1L]]),
  "Ora-F-Ada" %in% methods_for_scenario(4L, simulation4[[2L]]),
  !"Ora-F-Ada" %in% methods_for_scenario(4L, simulation4[[3L]]),
  !"Ora-F-Ada" %in% methods_for_scenario(4L, simulation4[[4L]])
)

cat(paste(
  "Validation tests passed: no calibration leakage, correct quantile logic,",
  "selection diagnostics, common candidate folds, disjoint TransLasso",
  "holdouts, an unlinked joint penalty grid, target-scaled lambda1,",
  "and valid oracle use.\n"
))
