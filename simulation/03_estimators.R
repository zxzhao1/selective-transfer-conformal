fit_source_initials <- function(sources, cfg, seed) {
  out <- vector("list", length(sources))
  for (k in seq_along(sources)) {
    source_k <- sources[[k]]
    fixed_lambda <- NULL
    if (!is.null(cfg$fixed_source_lambdas)) {
      fixed_lambda <- cfg$fixed_source_lambdas[k]
    }
    if (is.null(fixed_lambda) || is.na(fixed_lambda)) {
      grids <- make_lambda_grids(
        cfg,
        n_target = nrow(source_k$x),
        total_n = nrow(source_k$x),
        n_source = nrow(source_k$x)
      )
      fixed_lambda <- cv_lasso_lambda(
        source_k$x,
        source_k$y,
        grids$initial,
        cfg$inner_folds,
        seed + k,
        cfg$standardize
      )$lambda
    }
    out[[k]] <- fit_lasso_at_lambda(
      source_k$x, source_k$y, fixed_lambda, cfg$standardize
    )
  }
  out
}

get_fixed_method_parameters <- function(cfg, method) {
  parameters <- cfg$fixed_method_parameters
  if (is.null(parameters)) return(NULL)
  method_name <- as.character(method)
  selected <- parameters[
    parameters$method == method_name, , drop = FALSE
  ]
  if (nrow(selected) != 1L) {
    stop("Expected one fixed-parameter row for method ", method, ".",
         call. = FALSE)
  }
  selected
}

build_transfer_design <- function(target, sources) {
  p <- ncol(target$x)
  K <- length(sources)
  total_n <- nrow(target$x) + sum(vapply(sources, function(z) nrow(z$x), integer(1)))

  base_x <- matrix(0, total_n, p)
  source_indicators <- matrix(0, total_n, K)
  delta_x <- matrix(0, total_n, K * p)
  y <- numeric(total_n)

  cursor <- 0L
  target_rows <- seq_len(nrow(target$x))
  base_x[target_rows, ] <- target$x
  y[target_rows] <- target$y
  cursor <- nrow(target$x)

  for (k in seq_len(K)) {
    nk <- nrow(sources[[k]]$x)
    rows <- cursor + seq_len(nk)
    base_x[rows, ] <- sources[[k]]$x
    source_indicators[rows, k] <- 1
    block <- (k - 1L) * p + seq_len(p)
    delta_x[rows, block] <- sources[[k]]$x
    y[rows] <- sources[[k]]$y
    cursor <- cursor + nk
  }

  list(x = cbind(base_x, source_indicators, delta_x), y = y)
}

fit_transfer_fixed <- function(target, sources, cfg, source_initials,
                               lambda0, lambda1,
                               lambda_initial = NULL,
                               oracle_truth = NULL,
                               design = NULL,
                               target_initial = NULL) {
  p <- ncol(target$x)
  K <- length(sources)

  if (is.null(oracle_truth)) {
    if (is.null(target_initial)) {
      target_initial <- fit_lasso_at_lambda(
        target$x, target$y, lambda_initial, cfg$standardize
      )
    }
    w0 <- scad_weight(abs(target_initial$beta), lambda0, cfg$scad_a)
    wk <- lapply(seq_len(K), function(k) {
      initial_delta <- source_initials[[k]]$beta - target_initial$beta
      scad_weight(abs(initial_delta), lambda1, cfg$scad_a)
    })
  } else {
    tolerance <- sqrt(.Machine$double.eps)
    w0 <- as.numeric(abs(oracle_truth$beta0) <= tolerance)
    wk <- lapply(
      oracle_truth$deltas,
      function(delta) as.numeric(abs(delta) <= tolerance)
    )
  }

  if (is.null(design)) {
    design <- build_transfer_design(target, sources)
  }
  desired_penalty <- c(
    lambda0 * w0 / 2,
    rep(0, K),
    unlist(lapply(wk, function(w) lambda1 * w / 2), use.names = FALSE)
  )
  fit <- fit_weighted_lasso(
    design$x, design$y, desired_penalty, cfg$standardize
  )

  list(
    intercept = fit$intercept,
    beta = fit$coefficients[seq_len(p)],
    lambda0 = lambda0,
    lambda1 = lambda1,
    lambda_initial = lambda_initial,
    weights = list(target = w0, sources = wk),
    method = if (is.null(oracle_truth)) "F-AdaTrans" else "Ora-F-Ada"
  )
}

tune_transfer <- function(target, sources, cfg, source_initials, seed,
                          oracle_truth = NULL) {
  n_source <- mean(vapply(sources, function(z) nrow(z$x), integer(1)))
  total_n <- nrow(target$x) + sum(vapply(sources, function(z) nrow(z$x), integer(1)))
  grids <- make_lambda_grids(cfg, nrow(target$x), total_n, n_source)
  if (is.null(oracle_truth) && cfg$initial_tuning == "separate") {
    selected_initial <- cv_lasso_lambda(
      target$x,
      target$y,
      grids$initial,
      cfg$inner_folds,
      seed + 9000L,
      cfg$standardize
    )$lambda
    initial_grid <- selected_initial
  } else {
    initial_grid <- if (is.null(oracle_truth)) grids$initial else NA_real_
  }
  if (isTRUE(cfg$linked_transfer_penalties)) {
    penalty_pairs <- data.frame(
      lambda0 = grids$lambda0,
      lambda1 = grids$lambda1
    )
  } else {
    penalty_pairs <- expand.grid(
      lambda0 = grids$lambda0,
      lambda1 = grids$lambda1,
      KEEP.OUT.ATTRS = FALSE,
      stringsAsFactors = FALSE
    )
  }
  candidates <- merge(
    data.frame(lambda_initial = initial_grid),
    penalty_pairs,
    by = NULL
  )
  if (nrow(candidates) == 1L) {
    candidates$cv_error <- NA_real_
    return(candidates)
  }
  fold_id <- make_balanced_folds(
    nrow(target$x), cfg$inner_folds, seed, target$y
  )
  candidates$cv_error <- NA_real_

  folds <- sort(unique(fold_id))
  fold_cache <- lapply(folds, function(fold) {
    train <- fold_id != fold
    valid <- !train
    training_target <- list(
      x = target$x[train, , drop = FALSE],
      y = target$y[train]
    )
    initial_models <- NULL
    if (is.null(oracle_truth)) {
      initial_models <- lapply(initial_grid, function(lambda_initial) {
        fit_lasso_at_lambda(
          training_target$x,
          training_target$y,
          lambda_initial,
          cfg$standardize
        )
      })
    }
    list(
      train = train,
      valid = valid,
      target = training_target,
      design = build_transfer_design(training_target, sources),
      initial_models = initial_models
    )
  })

  for (g in seq_len(nrow(candidates))) {
    squared_errors <- numeric(nrow(target$x))
    initial_index <- if (is.null(oracle_truth)) {
      match(candidates$lambda_initial[g], initial_grid)
    } else {
      NA_integer_
    }
    for (fold_index in seq_along(folds)) {
      cache <- fold_cache[[fold_index]]
      model <- fit_transfer_fixed(
        target = cache$target,
        sources = sources,
        cfg = cfg,
        source_initials = source_initials,
        lambda0 = candidates$lambda0[g],
        lambda1 = candidates$lambda1[g],
        lambda_initial = candidates$lambda_initial[g],
        oracle_truth = oracle_truth,
        design = cache$design,
        target_initial = if (is.null(oracle_truth)) {
          cache$initial_models[[initial_index]]
        } else {
          NULL
        }
      )
      pred <- predict_linear_model(
        model, target$x[cache$valid, , drop = FALSE]
      )
      squared_errors[cache$valid] <- (target$y[cache$valid] - pred)^2
    }
    candidates$cv_error[g] <- mean(squared_errors)
  }

  candidates[which.min(candidates$cv_error), , drop = FALSE]
}

fit_tuned_transfer <- function(target, sources, cfg, source_initials, seed,
                               oracle_truth = NULL) {
  method <- if (is.null(oracle_truth)) "F-AdaTrans" else "Ora-F-Ada"
  fixed <- get_fixed_method_parameters(cfg, method)
  if (!is.null(fixed)) {
    model <- fit_transfer_fixed(
      target = target,
      sources = sources,
      cfg = cfg,
      source_initials = source_initials,
      lambda0 = fixed$lambda0,
      lambda1 = fixed$lambda1,
      lambda_initial = if (is.null(oracle_truth)) {
        fixed$lambda_initial
      } else {
        NA_real_
      },
      oracle_truth = oracle_truth
    )
    model$tuning <- data.frame(
      lambda_initial = model$lambda_initial,
      lambda0 = model$lambda0,
      lambda1 = model$lambda1,
      cv_error = NA_real_
    )
    return(model)
  }
  tuned <- tune_transfer(
    target, sources, cfg, source_initials, seed, oracle_truth
  )
  model <- fit_transfer_fixed(
    target = target,
    sources = sources,
    cfg = cfg,
    source_initials = source_initials,
    lambda0 = tuned$lambda0,
    lambda1 = tuned$lambda1,
    lambda_initial = tuned$lambda_initial,
    oracle_truth = oracle_truth
  )
  model$tuning <- tuned
  model
}

fit_tuned_target <- function(target, cfg, seed) {
  fixed <- get_fixed_method_parameters(cfg, "Target-only")
  if (!is.null(fixed)) {
    model <- fit_lasso_at_lambda(
      target$x, target$y, fixed$lambda_target, cfg$standardize
    )
    model$method <- "Target-only"
    model$tuning <- data.frame(
      lambda = fixed$lambda_target,
      cv_error = NA_real_
    )
    return(model)
  }
  grids <- make_lambda_grids(
    cfg, nrow(target$x), nrow(target$x), nrow(target$x)
  )
  tuned <- cv_lasso_lambda(
    target$x,
    target$y,
    grids$target,
    cfg$inner_folds,
    seed,
    cfg$standardize
  )
  model <- fit_lasso_at_lambda(
    target$x, target$y, tuned$lambda, cfg$standardize
  )
  model$method <- "Target-only"
  model$tuning <- tuned
  model
}

fit_pool_at_lambda <- function(target, sources, lambda, cfg) {
  x <- do.call(rbind, c(list(target$x), lapply(sources, `[[`, "x")))
  y <- unlist(c(list(target$y), lapply(sources, `[[`, "y")), use.names = FALSE)
  model <- fit_lasso_at_lambda(x, y, lambda, cfg$standardize)
  model$method <- "Complete pooling"
  model
}

fit_tuned_pool <- function(target, sources, cfg, seed) {
  fixed <- get_fixed_method_parameters(cfg, "Complete pooling")
  if (!is.null(fixed)) {
    model <- fit_pool_at_lambda(
      target, sources, fixed$lambda_target, cfg
    )
    model$tuning <- data.frame(
      lambda = fixed$lambda_target,
      cv_error = NA_real_
    )
    return(model)
  }
  total_n <- nrow(target$x) + sum(vapply(sources, function(z) nrow(z$x), integer(1)))
  grids <- make_lambda_grids(
    cfg,
    nrow(target$x),
    total_n,
    mean(vapply(sources, function(z) nrow(z$x), integer(1)))
  )
  fold_id <- make_balanced_folds(
    nrow(target$x), cfg$inner_folds, seed, target$y
  )
  cv_error <- numeric(length(grids$pool))

  for (g in seq_along(grids$pool)) {
    squared_errors <- numeric(nrow(target$x))
    for (fold in sort(unique(fold_id))) {
      train <- fold_id != fold
      valid <- !train
      model <- fit_pool_at_lambda(
        list(x = target$x[train, , drop = FALSE], y = target$y[train]),
        sources,
        grids$pool[g],
        cfg
      )
      pred <- predict_linear_model(model, target$x[valid, , drop = FALSE])
      squared_errors[valid] <- (target$y[valid] - pred)^2
    }
    cv_error[g] <- mean(squared_errors)
  }

  best <- which.min(cv_error)
  model <- fit_pool_at_lambda(target, sources, grids$pool[best], cfg)
  model$tuning <- data.frame(lambda = grids$pool, cv_error = cv_error)
  model
}
