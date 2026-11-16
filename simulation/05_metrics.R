evaluate_calibrated_method <- function(calibrated, target_test, beta0) {
  prediction <- predict_linear_model(calibrated$model, target_test$x)
  lower <- prediction - calibrated$q
  upper <- prediction + calibrated$q
  tuning <- calibrated$model$tuning
  lambda_initial <- NA_real_
  lambda0 <- NA_real_
  lambda1 <- NA_real_
  lambda_target <- NA_real_
  if (is.data.frame(tuning) && nrow(tuning) == 1L) {
    if ("lambda_initial" %in% names(tuning)) {
      lambda_initial <- tuning$lambda_initial[1L]
    }
    if ("lambda0" %in% names(tuning)) lambda0 <- tuning$lambda0[1L]
    if ("lambda1" %in% names(tuning)) lambda1 <- tuning$lambda1[1L]
  }
  if (!is.null(calibrated$model$lambda) &&
      length(calibrated$model$lambda) == 1L) {
    lambda_target <- calibrated$model$lambda
  }

  data.frame(
    method = calibrated$model$method,
    coverage = mean(target_test$y >= lower & target_test$y <= upper),
    interval_length = 2 * calibrated$q,
    prediction_mse = mean((target_test$y - prediction)^2),
    slope_error = sum((calibrated$model$beta - beta0)^2),
    lambda_initial = lambda_initial,
    lambda0 = lambda0,
    lambda1 = lambda1,
    lambda_target = lambda_target,
    trans_lambda_first = if (!is.null(calibrated$model$trans_lambda_first)) {
      calibrated$model$trans_lambda_first
    } else {
      NA_real_
    },
    trans_lambda_second = if (!is.null(calibrated$model$trans_lambda_second)) {
      calibrated$model$trans_lambda_second
    } else {
      NA_real_
    },
    selected_transfer = if (!is.null(calibrated$selected)) {
      as.integer(calibrated$selected == "transfer")
    } else {
      NA_integer_
    },
    selection_radius_transfer = if (!is.null(calibrated$r_transfer)) {
      calibrated$r_transfer
    } else {
      NA_real_
    },
    selection_radius_target = if (!is.null(calibrated$r_target)) {
      calibrated$r_target
    } else {
      NA_real_
    },
    selection_radius_margin = if (
      !is.null(calibrated$r_transfer) && !is.null(calibrated$r_target)
    ) {
      calibrated$r_target - calibrated$r_transfer
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
}

validate_simulation_data <- function(data, cfg) {
  stopifnot(
    ncol(data$target_fit$x) == cfg$p,
    nrow(data$target_fit$x) == data$setting$n_fit,
    nrow(data$target_cal$x) == data$setting$n_cal,
    nrow(data$target_test$x) == cfg$n_test,
    length(data$sources) == data$setting$K,
    length(data$truth$beta0) == cfg$p,
    all(vapply(data$sources, function(z) ncol(z$x), integer(1)) == cfg$p)
  )
  invisible(TRUE)
}
