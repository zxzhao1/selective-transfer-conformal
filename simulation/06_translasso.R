stable_translasso_aggregation <- function(
    B, X.test, y.test, total.step = 10L, selection = FALSE
) {
  B <- as.matrix(B)
  if (all(B == 0)) {
    theta <- numeric(ncol(B))
    theta[1L] <- 1
    return(list(theta = theta, beta = numeric(nrow(B)), beta.ew = numeric(nrow(B))))
  }
  softmax <- function(log_weight) {
    shifted <- log_weight - max(log_weight)
    weight <- exp(shifted)
    weight / sum(weight)
  }
  residual_sse <- colSums((as.numeric(y.test) - as.matrix(X.test) %*% B)^2)
  if (selection) {
    selected <- which.min(residual_sse)
    theta <- numeric(ncol(B))
    theta[selected] <- 1
    return(list(theta = theta, beta = B[, selected], beta.ew = NULL))
  }
  theta <- softmax(-residual_sse / 2)
  beta <- as.numeric(B %*% theta)
  beta_ew <- beta
  for (step in seq_len(total.step)) {
    disagreement <- colSums((
      as.numeric(as.matrix(X.test) %*% beta) - as.matrix(X.test) %*% B
    )^2)
    updated <- softmax(-residual_sse / 2 + disagreement / 8)
    beta <- as.numeric(B %*% updated / 4 + 3 * beta / 4)
    if (sum(abs(updated - theta)) < 1e-3) {
      theta <- updated
      break
    }
    theta <- updated
  }
  list(theta = theta, beta = beta, beta.ew = beta_ew)
}

load_translasso_author_code <- function(author_file) {
  if (is.null(author_file) || !file.exists(author_file)) {
    stop(
      "TransLasso author code was not found. Provide the path to ",
      "TransLasso-functions.R from https://github.com/saili0103/TransLasso.",
      call. = FALSE
    )
  }

  env <- new.env(parent = globalenv())
  env$glmnet <- glmnet::glmnet
  env$cv.glmnet <- glmnet::cv.glmnet
  sys.source(author_file, envir = env)
  if (!exists("Trans.lasso", envir = env, inherits = FALSE)) {
    stop("The supplied file does not define Trans.lasso().", call. = FALSE)
  }
  env$agg.fun <- stable_translasso_aggregation
  env
}

make_translasso_holdouts <- function(n_target, seed, y = NULL) {
  holdout_size <- floor(n_target / 3L)
  if (n_target - holdout_size < 8L || 2L * holdout_size > n_target) {
    stop("The target fitting sample is too small for the author implementation.",
         call. = FALSE)
  }
  set.seed(seed)
  if (!is.null(y)) {
    stopifnot(length(y) == n_target, !anyNA(y))
    if (length(unique(y)) <= max(20L, floor(n_target / 2L))) {
      partition <- make_stratified_partition(
        y, c(holdout_size, holdout_size, n_target - 2L * holdout_size)
      )
      return(list(
        first = which(partition == 1L),
        second = which(partition == 2L)
      ))
    }
  }
  target_order <- sample.int(n_target)
  list(
    first = target_order[seq_len(holdout_size)],
    second = target_order[holdout_size + seq_len(holdout_size)]
  )
}

select_translasso_lambda_constant <- function(x, y, n_vec, holdout) {
  x_train <- x[-holdout, , drop = FALSE]
  y_train <- y[-holdout]
  n_vec_train <- n_vec
  n_vec_train[1L] <- n_vec_train[1L] - length(holdout)
  target_rows <- seq_len(n_vec_train[1L])
  p <- ncol(x_train)
  lambda_grid <- seq(1, 0.1, length.out = 20L) *
    sqrt(2 * log(p) / n_vec_train[1L])
  tuned <- cv_lasso_lambda(
    x_train[target_rows, , drop = FALSE],
    y_train[target_rows],
    lambda_grid,
    min(8L, length(target_rows)),
    seed = 20261114L + sum(holdout),
    standardize = TRUE
  )
  tuned$lambda / sqrt(2 * log(p) / n_vec_train[1L])
}

fit_translasso_once_fixed <- function(author, x, y, n_vec, holdout,
                                      lambda_constant, l1 = TRUE) {
  x_holdout <- x[holdout, , drop = FALSE]
  y_holdout <- y[holdout]
  x <- x[-holdout, , drop = FALSE]
  y <- y[-holdout]
  n_vec[1L] <- n_vec[1L] - length(holdout)
  p <- ncol(x)
  target_rows <- author$ind.set(n_vec, 1L)
  source_count <- length(n_vec) - 1L

  source_scores <- numeric(source_count)
  for (k in seq_len(source_count)) {
    source_rows <- author$ind.set(n_vec, k + 1L)
    score_vector <-
      crossprod(x[source_rows, , drop = FALSE], y[source_rows]) / n_vec[k + 1L] -
      crossprod(x[target_rows, , drop = FALSE], y[target_rows]) / n_vec[1L]
    margin <- sort(abs(score_vector), decreasing = TRUE)[
      seq_len(round(n_vec[1L] / 3L))
    ]
    source_scores[k] <- sum(margin^2)
  }
  source_ranks <- rank(source_scores)
  rank_levels <- unique(source_ranks)
  transfer_sets <- lapply(
    rank_levels,
    function(level) which(source_ranks <= level)
  )
  transfer_sets <- unique(transfer_sets)

  target_lambda <- lambda_constant * sqrt(2 * log(p) / n_vec[1L])
  target_lambda_path <- seq(1, 0.1, length.out = 20L) *
    sqrt(2 * log(p) / n_vec[1L])
  target_fit <- glmnet::glmnet(
    x[target_rows, , drop = FALSE],
    y[target_rows],
    lambda = target_lambda_path
  )
  target_beta <- as.numeric(predict(
    target_fit, s = target_lambda, type = "coefficients"
  )[-1L])
  candidates <- list(target_beta)
  for (k in seq_along(transfer_sets)) {
    candidates[[k + 1L]] <- author$las.kA(
      X = x,
      y = y,
      A0 = transfer_sets[[k]],
      n.vec = n_vec,
      lam.const = lambda_constant,
      l1 = l1
    )$beta.kA
  }
  candidates <- candidates[!duplicated(candidates)]
  candidate_matrix <- as.matrix(as.data.frame(candidates))
  aggregation <- author$agg.fun(
    B = candidate_matrix, X.test = x_holdout, y.test = y_holdout
  )
  list(
    beta.hat = aggregation$beta,
    theta.hat = aggregation$theta,
    rank.pi = source_ranks,
    lam.const = lambda_constant
  )
}

fit_translasso_author <- function(data, cfg, author_file, seed) {
  author <- load_translasso_author_code(author_file)
  x <- do.call(
    rbind,
    c(list(data$target_fit$x), lapply(data$sources, `[[`, "x"))
  )
  y <- unlist(
    c(list(data$target_fit$y), lapply(data$sources, `[[`, "y")),
    use.names = FALSE
  )
  n_vec <- c(
    nrow(data$target_fit$x),
    vapply(data$sources, function(z) nrow(z$x), integer(1))
  )

  holdouts <- make_translasso_holdouts(n_vec[1L], seed, data$target_fit$y)
  fixed <- get_fixed_method_parameters(cfg, "TransLasso")
  if (is.null(fixed)) {
    lambda_first <- select_translasso_lambda_constant(
      x, y, n_vec, holdouts$first
    )
    lambda_second <- select_translasso_lambda_constant(
      x, y, n_vec, holdouts$second
    )
    first <- fit_translasso_once_fixed(
      author, x, y, n_vec, holdouts$first, lambda_first, l1 = TRUE
    )
    second <- fit_translasso_once_fixed(
      author, x, y, n_vec, holdouts$second, lambda_second, l1 = TRUE
    )
  } else {
    first <- fit_translasso_once_fixed(
      author, x, y, n_vec, holdouts$first,
      fixed$trans_lambda_first, l1 = TRUE
    )
    second <- fit_translasso_once_fixed(
      author, x, y, n_vec, holdouts$second,
      fixed$trans_lambda_second, l1 = TRUE
    )
  }

  list(
    intercept = 0,
    beta = as.numeric((first$beta.hat + second$beta.hat) / 2),
    method = "TransLasso",
    author_file = normalizePath(author_file),
    author_file_md5 = unname(tools::md5sum(author_file)),
    holdouts = holdouts,
    trans_lambda_first = if (!is.null(first$lam.const)) {
      first$lam.const
    } else {
      NA_real_
    },
    trans_lambda_second = if (!is.null(second$lam.const)) {
      second$lam.const
    } else {
      NA_real_
    }
  )
}

fit_and_calibrate_translasso <- function(data, cfg, author_file, seed) {
  model <- fit_translasso_author(data, cfg, author_file, seed)
  calibrate_model(model, data$target_cal, cfg$alpha)
}
