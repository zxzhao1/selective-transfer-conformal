make_scenarios <- function(simulation, cfg) {
  simulation <- as.integer(simulation)

  if (simulation == 1L) {
    values <- c(0, 0.3, 0.6, 0.9)
    return(lapply(seq_along(values), function(j) list(
      id = sprintf("sim1_h_%s", values[j]),
      label = sprintf("h_delta=%s", values[j]),
      h_delta = values[j],
      s_delta = min(25L, cfg$p),
      K = cfg$K,
      n_source = cfg$n_source,
      n_fit = cfg$n_fit,
      error_case = "gaussian"
    )))
  }

  if (simulation == 2L) {
    values <- unique(pmin(c(0L, 10L, 25L, 50L, 100L), cfg$p))
    return(lapply(seq_along(values), function(j) list(
      id = sprintf("sim2_s_%s", values[j]),
      label = sprintf("s_delta=%s", values[j]),
      h_delta = 0.6,
      s_delta = values[j],
      K = cfg$K,
      n_source = cfg$n_source,
      n_fit = cfg$n_fit,
      error_case = "gaussian"
    )))
  }

  if (simulation == 3L) {
    source_sizes <- if (cfg$mode == "full") c(100L, 250L, 500L) else c(30L, 50L, 80L)
    source_counts <- if (cfg$mode == "full") c(2L, 4L, 8L) else c(2L, 3L, 4L)
    target_sizes <- if (cfg$mode == "full") c(30L, 50L, 100L) else c(24L, 30L, 42L)

    scenarios <- list()
    for (value in source_sizes) {
      scenarios[[length(scenarios) + 1L]] <- list(
        id = sprintf("sim3_nsource_%s", value),
        label = sprintf("n_source=%s", value),
        h_delta = cfg$h_delta,
        s_delta = cfg$s_delta,
        K = cfg$K,
        n_source = value,
        n_fit = cfg$n_fit,
        error_case = "gaussian"
      )
    }
    for (value in source_counts) {
      scenarios[[length(scenarios) + 1L]] <- list(
        id = sprintf("sim3_K_%s", value),
        label = sprintf("K=%s", value),
        h_delta = cfg$h_delta,
        s_delta = cfg$s_delta,
        K = value,
        n_source = cfg$n_source,
        n_fit = cfg$n_fit,
        error_case = "gaussian"
      )
    }
    for (value in target_sizes) {
      scenarios[[length(scenarios) + 1L]] <- list(
        id = sprintf("sim3_nfit_%s", value),
        label = sprintf("n_fit=%s", value),
        h_delta = cfg$h_delta,
        s_delta = cfg$s_delta,
        K = cfg$K,
        n_source = cfg$n_source,
        n_fit = value,
        error_case = "gaussian"
      )
    }
    return(scenarios)
  }

  if (simulation == 4L) {
    return(list(
      list(
        id = "sim4_complete",
        label = "complete transferability",
        h_delta = 0,
        s_delta = 0L,
        K = cfg$K,
        n_source = cfg$n_source,
        n_fit = cfg$n_fit,
        error_case = "gaussian"
      ),
      list(
        id = "sim4_sparse",
        label = "sparse feature-wise heterogeneity",
        h_delta = 0.6,
        s_delta = min(25L, cfg$p),
        K = cfg$K,
        n_source = cfg$n_source,
        n_fit = cfg$n_fit,
        error_case = "gaussian"
      ),
      list(
        id = "sim4_heterogeneous",
        label = "heterogeneous source transferability",
        h_delta = rep_len(c(0, 0.3, 0.6, 0.9), cfg$K),
        s_delta = rep_len(pmin(c(0L, 25L, 100L, 500L), cfg$p), cfg$K),
        K = cfg$K,
        n_source = cfg$n_source,
        n_fit = cfg$n_fit,
        error_case = "gaussian"
      ),
      list(
        id = "sim4_dense",
        label = "dense feature-wise heterogeneity",
        h_delta = 0.9,
        s_delta = cfg$p,
        K = cfg$K,
        n_source = cfg$n_source,
        n_fit = cfg$n_fit,
        error_case = "gaussian"
      )
    ))
  }

  if (simulation == 5L) {
    cases <- c("gaussian", "variance_shift", "heavy_tail")
    return(lapply(cases, function(error_case) list(
      id = paste0("sim5_", error_case),
      label = error_case,
      h_delta = cfg$h_delta,
      s_delta = cfg$s_delta,
      K = cfg$K,
      n_source = cfg$n_source,
      n_fit = cfg$n_fit,
      error_case = error_case
    )))
  }

  stop("simulation must be an integer from 1 to 5.", call. = FALSE)
}

methods_for_simulation <- function(simulation) {
  if (as.integer(simulation) <= 3L) {
    c("Proposed", "F-AdaTrans", "Target-only")
  } else {
    c(
      "Proposed", "F-AdaTrans", "Target-only", "Complete pooling",
      "TransLasso", "Ora-F-Ada"
    )
  }
}

oracle_is_well_defined <- function(scenario) {
  K <- as.integer(scenario$K)
  s_delta <- as.integer(rep_len(scenario$s_delta, K))
  n_source <- as.integer(rep_len(scenario$n_source, K))
  all(s_delta < n_source)
}

methods_for_scenario <- function(simulation, scenario) {
  methods <- methods_for_simulation(simulation)
  if ("Ora-F-Ada" %in% methods && !oracle_is_well_defined(scenario)) {
    methods <- setdiff(methods, "Ora-F-Ada")
  }
  methods
}
