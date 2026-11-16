script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1L])))

for (file in c(
  "00_config.R", "01_utils.R", "02_dgp.R", "03_estimators.R",
  "04_proposed.R", "05_metrics.R", "06_translasso.R", "07_scenarios.R"
)) {
  source(file.path(script_dir, file), local = FALSE)
}

require_simulation_packages()
args <- commandArgs(trailingOnly = TRUE)
pilot_file <- if (length(args) >= 1L) args[1L] else file.path(
  script_dir, "results", "simulation_4_full_compact_pilot10.csv"
)
output_file <- if (length(args) >= 2L) args[2L] else file.path(
  script_dir, "results", "simulation_4_pilot10_fixed_parameters.csv"
)
simulation_number <- if (length(args) >= 3L) as.integer(args[3L]) else 4L

pilot <- data.table::fread(pilot_file)
pilot_keys <- unique(pilot[, .(scenario_number, replication, seed)])
expected_scenarios <- length(make_scenarios(simulation_number, make_config("full")))
if (any(pilot_keys$replication > 10L) ||
    any(pilot_keys$replication < 1L) ||
    nrow(pilot_keys) != 10L * expected_scenarios) {
  stop("The pilot file must contain replications 1--10 in every scenario.",
       call. = FALSE)
}

cfg <- make_config("full")
cfg$initial_tuning <- "joint"
cfg$initial_lambda_multipliers <- c(0.5, 1)
cfg$target_lambda_multipliers <- c(0.5, 1)
cfg$pool_lambda_multipliers <- c(0.5, 1, 2, 4, 8)
cfg$transfer_lambda_multipliers <- c(0.5, 1, 2, 4, 8)
cfg$linked_transfer_penalties <- TRUE
scenarios <- make_scenarios(simulation_number, cfg)

source_records <- list()
translasso_records <- list()
record_index <- 0L
for (scenario_index in seq_along(scenarios)) {
  scenario <- scenarios[[scenario_index]]
  local_cfg <- cfg
  local_cfg$K <- as.integer(scenario$K)
  local_cfg$scad_a <- 3.7 * local_cfg$K / 2
  local_cfg$n_fit <- as.integer(scenario$n_fit)
  local_cfg$n_source <- as.integer(scenario$n_source)
  local_cfg <- validate_config(local_cfg)

  for (replication_index in seq_len(10L)) {
    seed <- pilot_keys[
      scenario_number == scenario_index & replication == replication_index,
      seed
    ]
    if (length(seed) != 1L) {
      stop("Could not resolve a unique pilot seed.", call. = FALSE)
    }
    data <- generate_simulation_data(
      local_cfg,
      seed = seed,
      h_delta = scenario$h_delta,
      s_delta = scenario$s_delta,
      K = scenario$K,
      n_source = scenario$n_source,
      n_fit = scenario$n_fit,
      n_cal = local_cfg$n_cal,
      error_case = scenario$error_case
    )
    source_initials <- fit_source_initials(
      data$sources, local_cfg, seed + 10L
    )
    for (source_index in seq_along(source_initials)) {
      record_index <- record_index + 1L
      source_records[[record_index]] <- data.table::data.table(
        scenario_number = scenario_index,
        replication = replication_index,
        source_index = source_index,
        lambda_initial = source_initials[[source_index]]$lambda
      )
    }

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
    holdouts <- make_translasso_holdouts(n_vec[1L], seed + 500L)
    lambda_first <- select_translasso_lambda_constant(
      x, y, n_vec, holdouts$first
    )
    lambda_second <- select_translasso_lambda_constant(
      x, y, n_vec, holdouts$second
    )
    translasso_records[[length(translasso_records) + 1L]] <-
      data.table::data.table(
        scenario_number = scenario_index,
        replication = replication_index,
        trans_lambda_first = lambda_first,
        trans_lambda_second = lambda_second
      )
  }
}

source_means <- data.table::rbindlist(source_records)[
  , .(lambda_initial = mean(lambda_initial)),
  by = .(scenario_number, source_index)
]
translasso_means <- data.table::rbindlist(translasso_records)[
  , .(
    trans_lambda_first = mean(trans_lambda_first),
    trans_lambda_second = mean(trans_lambda_second)
  ),
  by = scenario_number
]

numeric_mean <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}
method_means <- pilot[method %in% c(
  "F-AdaTrans", "Target-only", "Complete pooling", "Ora-F-Ada"
), .(
  lambda_initial = numeric_mean(lambda_initial),
  lambda0 = numeric_mean(lambda0),
  lambda1 = numeric_mean(lambda1),
  lambda_target = numeric_mean(lambda_target)
), by = .(scenario_number, method)]

translasso_rows <- translasso_means[, .(
  scenario_number,
  method = "TransLasso",
  lambda_initial = NA_real_,
  lambda0 = NA_real_,
  lambda1 = NA_real_,
  lambda_target = NA_real_,
  trans_lambda_first,
  trans_lambda_second
)]
method_means[, `:=`(
  trans_lambda_first = NA_real_,
  trans_lambda_second = NA_real_
)]
method_rows <- data.table::rbindlist(
  list(method_means, translasso_rows), use.names = TRUE, fill = TRUE
)
method_rows[, `:=`(
  parameter_role = "method",
  source_index = NA_integer_
)]

source_rows <- source_means[, .(
  scenario_number,
  parameter_role = "source_initial",
  method = "Source-initial",
  source_index,
  lambda_initial,
  lambda0 = NA_real_,
  lambda1 = NA_real_,
  lambda_target = NA_real_,
  trans_lambda_first = NA_real_,
  trans_lambda_second = NA_real_
)]

output <- data.table::rbindlist(
  list(method_rows, source_rows), use.names = TRUE, fill = TRUE
)
data.table::setcolorder(output, c(
  "scenario_number", "parameter_role", "method", "source_index",
  "lambda_initial", "lambda0", "lambda1", "lambda_target",
  "trans_lambda_first", "trans_lambda_second"
))
data.table::setorder(output, scenario_number, parameter_role, method, source_index)
data.table::fwrite(output, output_file)
print(output)
message("Saved fixed parameters to ", normalizePath(output_file))
