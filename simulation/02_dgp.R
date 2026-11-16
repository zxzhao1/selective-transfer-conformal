draw_covariates <- function(n, sigma) {
  z <- matrix(stats::rnorm(n * ncol(sigma)), nrow = n)
  z %*% chol(sigma)
}

draw_source_errors <- function(n, source_index, error_case) {
  if (error_case == "gaussian") {
    return(stats::rnorm(n))
  }
  if (error_case == "variance_shift") {
    source_sd <- c(0.5, 1, 1.5, 2)
    return(stats::rnorm(n, sd = source_sd[(source_index - 1L) %% 4L + 1L]))
  }
  if (error_case == "heavy_tail") {
    return(sqrt(3 / 5) * stats::rt(n, df = 5))
  }
  stop("Unknown error case: ", error_case, call. = FALSE)
}

make_contrast_support <- function(p, s_delta, source_index) {
  if (s_delta == 0L) {
    return(integer(0))
  }
  if (source_index %% 2L == 1L) {
    seq_len(s_delta)
  } else {
    seq.int(p - s_delta + 1L, p)
  }
}

generate_simulation_data <- function(cfg, seed = cfg$seed,
                                     h_delta = cfg$h_delta,
                                     s_delta = cfg$s_delta,
                                     K = cfg$K,
                                     n_source = cfg$n_source,
                                     n_fit = cfg$n_fit,
                                     n_cal = cfg$n_cal,
                                     error_case = "gaussian") {
  set.seed(seed)
  p <- cfg$p
  h_by_source <- rep_len(h_delta, K)
  s_by_source <- as.integer(rep_len(s_delta, K))
  n_by_source <- as.integer(rep_len(n_source, K))

  if (any(s_by_source < 0L | s_by_source > p)) {
    stop("Every s_delta must lie between 0 and p.", call. = FALSE)
  }

  beta0 <- numeric(p)
  beta0[seq_len(cfg$s0)] <- cfg$beta_signal

  n_target <- n_fit + n_cal
  x_target <- matrix(stats::rnorm(n_target * p), nrow = n_target)
  y_target <- drop(x_target %*% beta0 + stats::rnorm(n_target))
  split_order <- sample.int(n_target)
  fit_index <- split_order[seq_len(n_fit)]
  cal_index <- split_order[n_fit + seq_len(n_cal)]

  target_fit <- list(
    x = x_target[fit_index, , drop = FALSE],
    y = y_target[fit_index]
  )
  target_cal <- list(
    x = x_target[cal_index, , drop = FALSE],
    y = y_target[cal_index]
  )

  x_test <- matrix(stats::rnorm(cfg$n_test * p), nrow = cfg$n_test)
  y_test <- drop(x_test %*% beta0 + stats::rnorm(cfg$n_test))
  target_test <- list(x = x_test, y = y_test)

  sources <- vector("list", K)
  deltas <- vector("list", K)
  sigmas <- vector("list", K)

  for (k in seq_len(K)) {
    a <- matrix(
      ifelse(stats::runif(p * p) < 0.3, 0.3, 0),
      nrow = p,
      ncol = p
    )
    sigma_k <- crossprod(a) + diag(p)
    x_k <- draw_covariates(n_by_source[k], sigma_k)
    support_k <- make_contrast_support(p, s_by_source[k], k)
    delta_k <- numeric(p)
    if (length(support_k) > 0L && h_by_source[k] > 0) {
      delta_k[support_k] <- stats::rnorm(
        length(support_k),
        mean = h_by_source[k],
        sd = h_by_source[k] / 3
      )
    }
    errors_k <- draw_source_errors(n_by_source[k], k, error_case)
    y_k <- drop(x_k %*% (beta0 + delta_k) + errors_k)

    sources[[k]] <- list(x = x_k, y = y_k)
    deltas[[k]] <- delta_k
    sigmas[[k]] <- sigma_k
  }

  list(
    target_fit = target_fit,
    target_cal = target_cal,
    target_test = target_test,
    sources = sources,
    truth = list(beta0 = beta0, deltas = deltas, sigmas = sigmas),
    setting = list(
      h_delta = h_by_source,
      s_delta = s_by_source,
      K = K,
      n_source = n_by_source,
      n_fit = n_fit,
      n_cal = n_cal,
      error_case = error_case
    )
  )
}
