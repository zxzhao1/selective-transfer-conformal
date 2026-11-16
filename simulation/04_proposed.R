fit_proposed <- function(target_fit, target_cal, sources, cfg,
                         source_initials, seed) {
  fold_id <- make_balanced_folds(
    nrow(target_fit$x), cfg$outer_folds, seed, target_fit$y
  )
  transfer_errors <- numeric(nrow(target_fit$x))
  target_errors <- numeric(nrow(target_fit$x))

  for (fold in sort(unique(fold_id))) {
    train <- fold_id != fold
    valid <- !train
    outer_train <- list(
      x = target_fit$x[train, , drop = FALSE],
      y = target_fit$y[train]
    )

    transfer_model <- fit_tuned_transfer(
      outer_train,
      sources,
      cfg,
      source_initials,
      seed + 100L + fold
    )
    target_model <- fit_tuned_target(
      outer_train,
      cfg,
      seed + 200L + fold
    )

    transfer_pred <- predict_linear_model(
      transfer_model, target_fit$x[valid, , drop = FALSE]
    )
    target_pred <- predict_linear_model(
      target_model, target_fit$x[valid, , drop = FALSE]
    )
    transfer_errors[valid] <- abs(target_fit$y[valid] - transfer_pred)
    target_errors[valid] <- abs(target_fit$y[valid] - target_pred)
  }

  r_transfer <- finite_sample_quantile(transfer_errors, cfg$alpha)
  r_target <- finite_sample_quantile(target_errors, cfg$alpha)
  selected <- if (r_transfer < r_target) "transfer" else "target"

  if (selected == "transfer") {
    final_model <- fit_tuned_transfer(
      target_fit,
      sources,
      cfg,
      source_initials,
      seed + 300L
    )
  } else {
    final_model <- fit_tuned_target(target_fit, cfg, seed + 300L)
  }
  final_model$method <- "Proposed"
  calibrated <- calibrate_model(final_model, target_cal, cfg$alpha)
  calibrated$selected <- selected
  calibrated$r_transfer <- r_transfer
  calibrated$r_target <- r_target
  calibrated$outer_errors <- list(
    transfer = transfer_errors,
    target = target_errors
  )
  calibrated
}

fit_and_calibrate <- function(method, data, cfg, source_initials, seed) {
  model <- switch(
    method,
    "F-AdaTrans" = fit_tuned_transfer(
      data$target_fit, data$sources, cfg, source_initials, seed
    ),
    "Target-only" = fit_tuned_target(data$target_fit, cfg, seed),
    "Complete pooling" = fit_tuned_pool(
      data$target_fit, data$sources, cfg, seed
    ),
    "Ora-F-Ada" = fit_tuned_transfer(
      data$target_fit,
      data$sources,
      cfg,
      source_initials,
      seed,
      oracle_truth = data$truth
    ),
    stop("Unknown method: ", method, call. = FALSE)
  )
  calibrate_model(model, data$target_cal, cfg$alpha)
}
