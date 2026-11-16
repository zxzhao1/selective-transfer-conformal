require_simulation_packages <- function() {
  required <- c("glmnet", "Matrix", "data.table", "ggplot2", "scales")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop("Missing R packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

make_stratified_partition <- function(strata, capacities, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  capacities <- as.integer(capacities)
  stopifnot(
    length(strata) == sum(capacities), all(capacities >= 0L), !anyNA(strata)
  )
  groups <- lapply(
    split(seq_along(strata), as.character(strata)),
    function(rows) rows[sample.int(length(rows))]
  )
  group_order <- order(lengths(groups), runif(length(groups)))
  groups <- groups[group_order]
  proportions <- capacities / sum(capacities)
  assigned <- integer(length(capacities))
  partition <- integer(length(strata))
  processed <- 0L

  for (rows in groups) {
    for (row in rows) {
      feasible <- which(assigned < capacities)
      deficit <- proportions * (processed + 1L) - assigned
      best <- feasible[deficit[feasible] == max(deficit[feasible])]
      chosen <- best[sample.int(length(best), 1L)]
      partition[row] <- chosen
      assigned[chosen] <- assigned[chosen] + 1L
      processed <- processed + 1L
    }
  }
  stopifnot(identical(as.integer(tabulate(partition)), capacities))
  partition
}

make_balanced_folds <- function(n, nfolds, seed, strata = NULL) {
  nfolds <- min(as.integer(nfolds), as.integer(n))
  if (nfolds < 2L) {
    stop("At least two folds are required.", call. = FALSE)
  }
  set.seed(seed)
  if (!is.null(strata)) {
    stopifnot(length(strata) == n, !anyNA(strata))
    if (length(unique(strata)) <= max(20L, floor(n / 2L))) {
      capacities <- tabulate(
        rep(seq_len(nfolds), length.out = n), nbins = nfolds
      )
      return(make_stratified_partition(strata, capacities))
    }
  }
  sample(rep(seq_len(nfolds), length.out = n))
}

fit_lasso_at_lambda <- function(x, y, lambda, standardize = TRUE) {
  y <- as.numeric(y)
  if (length(unique(y[is.finite(y)])) < 2L || stats::sd(y) <= 1e-12) {
    return(list(
      intercept = mean(y),
      beta = numeric(ncol(x)),
      lambda = lambda,
      fit = NULL,
      method = "lasso"
    ))
  }
  fit <- glmnet::glmnet(
    x = as.matrix(x),
    y = y,
    alpha = 1,
    lambda = lambda,
    intercept = TRUE,
    standardize = standardize
  )
  coefs <- as.numeric(stats::coef(fit, s = lambda))
  list(
    intercept = coefs[1L],
    beta = coefs[-1L],
    lambda = lambda,
    fit = fit,
    method = "lasso"
  )
}

predict_linear_model <- function(model, newx) {
  drop(model$intercept + as.matrix(newx) %*% model$beta)
}

cv_lasso_lambda <- function(x, y, lambdas, nfolds, seed,
                            standardize = TRUE) {
  if (length(lambdas) == 1L) {
    return(list(
      lambda = lambdas[1L],
      cv_error = NA_real_,
      table = data.frame(lambda = lambdas[1L], cv_error = NA_real_)
    ))
  }
  fold_id <- make_balanced_folds(nrow(x), nfolds, seed, y)
  cv_error <- numeric(length(lambdas))

  for (g in seq_along(lambdas)) {
    squared_errors <- numeric(nrow(x))
    for (fold in sort(unique(fold_id))) {
      train <- fold_id != fold
      valid <- !train
      model <- fit_lasso_at_lambda(
        x[train, , drop = FALSE], y[train], lambdas[g], standardize
      )
      pred <- predict_linear_model(model, x[valid, , drop = FALSE])
      squared_errors[valid] <- (y[valid] - pred)^2
    }
    cv_error[g] <- mean(squared_errors)
  }

  best <- which.min(cv_error)
  list(
    lambda = lambdas[best],
    cv_error = cv_error[best],
    table = data.frame(lambda = lambdas, cv_error = cv_error)
  )
}

scad_weight <- function(abs_value, lambda, a) {
  if (lambda <= 0 || a <= 2) {
    stop("SCAD requires lambda > 0 and a > 2.", call. = FALSE)
  }
  out <- numeric(length(abs_value))
  first <- abs_value <= lambda
  middle <- abs_value > lambda & abs_value <= a * lambda
  out[first] <- 1
  out[middle] <- (a * lambda - abs_value[middle]) / ((a - 1) * lambda)
  out
}

fit_weighted_lasso <- function(x, y, desired_penalty,
                               standardize = TRUE) {
  desired_penalty <- as.numeric(desired_penalty)
  if (length(desired_penalty) != ncol(x)) {
    stop("Penalty vector and design matrix have incompatible dimensions.",
         call. = FALSE)
  }
  if (any(desired_penalty < 0) || sum(desired_penalty) <= 0) {
    stop("The weighted Lasso requires a nonzero, nonnegative penalty.",
         call. = FALSE)
  }

  # glmnet rescales penalty.factor to sum to the number of columns.
  # This lambda value makes lambda * penalty.factor equal to the desired
  # coefficient-specific penalty in glmnet's 1/(2N) loss convention.
  lambda_glmnet <- sum(desired_penalty) / length(desired_penalty)
  fit <- glmnet::glmnet(
    x = as.matrix(x),
    y = as.numeric(y),
    alpha = 1,
    lambda = lambda_glmnet,
    penalty.factor = desired_penalty,
    intercept = TRUE,
    standardize = standardize
  )
  coefs <- as.numeric(stats::coef(fit, s = lambda_glmnet))
  list(intercept = coefs[1L], coefficients = coefs[-1L], fit = fit)
}

finite_sample_quantile <- function(scores, alpha) {
  scores <- sort(as.numeric(scores))
  rank <- ceiling((length(scores) + 1) * (1 - alpha))
  if (rank > length(scores)) Inf else scores[rank]
}

calibrate_model <- function(model, calibration, alpha) {
  residuals <- abs(
    calibration$y - predict_linear_model(model, calibration$x)
  )
  list(
    model = model,
    q = finite_sample_quantile(residuals, alpha),
    residuals = residuals
  )
}
