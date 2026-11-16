script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1L])))

for (file in c(
  "00_config.R",
  "01_utils.R",
  "02_dgp.R",
  "03_estimators.R",
  "04_proposed.R",
  "05_metrics.R"
)) {
  source(file.path(script_dir, file), local = FALSE)
}

require_simulation_packages()
cfg <- make_config("smoke")
data <- generate_simulation_data(cfg, seed = cfg$seed)
validate_simulation_data(data, cfg)

source_initials <- fit_source_initials(data$sources, cfg, cfg$seed + 10L)
proposed <- fit_proposed(
  data$target_fit,
  data$target_cal,
  data$sources,
  cfg,
  source_initials,
  cfg$seed + 20L
)

methods <- c("F-AdaTrans", "Target-only", "Complete pooling", "Ora-F-Ada")
fitted <- lapply(seq_along(methods), function(j) {
  fit_and_calibrate(
    methods[j], data, cfg, source_initials, cfg$seed + 1000L + j
  )
})
names(fitted) <- methods
fitted <- c(list(Proposed = proposed), fitted)

results <- data.table::rbindlist(lapply(fitted, function(object) {
  evaluate_calibrated_method(object, data$target_test, data$truth$beta0)
}))

stopifnot(
  nrow(results) == 5L,
  all(is.finite(results$coverage)),
  all(results$coverage >= 0 & results$coverage <= 1),
  all(results$interval_length > 0),
  all(is.finite(results$prediction_mse)),
  all(is.finite(results$slope_error)),
  proposed$selected %in% c("transfer", "target"),
  length(proposed$outer_errors$transfer) == cfg$n_fit,
  length(proposed$residuals) == cfg$n_cal
)

print(results)
cat("Smoke test passed. Selected estimator:", proposed$selected, "\n")
