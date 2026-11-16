make_config <- function(mode = c("smoke", "full")) {
  mode <- match.arg(mode)

  if (mode == "smoke") {
    cfg <- list(
      mode = mode,
      p = 30L,
      K = 2L,
      n_fit = 30L,
      n_cal = 40L,
      n_source = 50L,
      n_test = 100L,
      s0 = 5L,
      beta_signal = 0.3,
      h_delta = 0.6,
      s_delta = 6L,
      alpha = 0.1,
      outer_folds = 3L,
      inner_folds = 3L,
      replications = 2L,
      lambda_multipliers = c(0.7, 1.4),
      linked_transfer_penalties = FALSE,
      initial_tuning = "joint",
      standardize = TRUE,
      seed = 20261113L
    )
  } else {
    cfg <- list(
      mode = mode,
      p = 500L,
      K = 4L,
      n_fit = 50L,
      n_cal = 100L,
      n_source = 250L,
      n_test = 2000L,
      s0 = 8L,
      beta_signal = 0.3,
      h_delta = 0.6,
      s_delta = 25L,
      alpha = 0.1,
      outer_folds = 5L,
      inner_folds = 3L,
      replications = 200L,
      lambda_multipliers = c(0.5, 1, 2),
      linked_transfer_penalties = FALSE,
      initial_tuning = "joint",
      standardize = TRUE,
      seed = 20261113L
    )
  }

  cfg$implementation_version <- "2026-11-13-v6-fixed-penalties"
  cfg$scad_a <- 3.7 * cfg$K / 2
  validate_config(cfg)
}

validate_config <- function(cfg) {
  stopifnot(
    cfg$p >= 2L,
    cfg$K >= 1L,
    cfg$n_fit >= cfg$outer_folds,
    cfg$n_fit >= cfg$inner_folds,
    cfg$n_cal >= 1L,
    cfg$n_source >= cfg$inner_folds,
    cfg$n_test >= 1L,
    cfg$s0 <= cfg$p,
    cfg$s_delta <= cfg$p,
    cfg$alpha > 0,
    cfg$alpha < 1,
    cfg$scad_a > 2,
    is.character(cfg$implementation_version),
    length(cfg$implementation_version) == 1L,
    nzchar(cfg$implementation_version),
    all(cfg$lambda_multipliers > 0),
    cfg$initial_tuning %in% c("joint", "separate")
  )
  cfg
}

make_lambda_grids <- function(cfg, n_target, total_n, n_source) {
  p_log <- log(max(cfg$p, 2L))
  if (!is.null(cfg$fixed_penalty_multipliers)) {
    mult_initial <- cfg$fixed_penalty_multipliers$initial
    mult_target <- cfg$fixed_penalty_multipliers$target
    mult_pool <- cfg$fixed_penalty_multipliers$pool
    mult_lambda0 <- cfg$fixed_penalty_multipliers$lambda0
    mult_lambda1 <- cfg$fixed_penalty_multipliers$lambda1
  } else {
    mult_initial <- if (is.null(cfg$initial_lambda_multipliers)) {
      cfg$lambda_multipliers
    } else {
      cfg$initial_lambda_multipliers
    }
    mult_target <- if (is.null(cfg$target_lambda_multipliers)) {
      cfg$lambda_multipliers
    } else {
      cfg$target_lambda_multipliers
    }
    mult_pool <- if (is.null(cfg$pool_lambda_multipliers)) {
      cfg$lambda_multipliers
    } else {
      cfg$pool_lambda_multipliers
    }
    transfer_multipliers <- if (is.null(cfg$transfer_lambda_multipliers)) {
      cfg$lambda_multipliers
    } else {
      cfg$transfer_lambda_multipliers
    }
    mult_lambda0 <- mult_lambda1 <- transfer_multipliers
  }

  list(
    initial = mult_initial * sqrt(p_log / n_target),
    target = mult_target * sqrt(p_log / n_target),
    pool = mult_pool * sqrt(p_log / total_n),
    lambda0 = mult_lambda0 * sqrt(p_log / total_n),
    lambda1 = mult_lambda1 * (n_source / total_n) *
      sqrt(p_log / n_target)
  )
}
