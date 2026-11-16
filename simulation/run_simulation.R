parse_cli <- function(args) {
  defaults <- list(
    simulation = "1",
    mode = "smoke",
    replications = NA_character_,
    scenario = "all",
    translasso = NA_character_,
    overwrite = "false",
    workers = "1",
    run_id = "",
    resume = "true",
    initial_tuning = "joint",
    grid = "standard",
    seed_offset = "0",
    fixed_params = NA_character_
  )
  for (arg in args) {
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1L]]
    if (length(pieces) != 2L || !pieces[1L] %in% names(defaults)) {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
    defaults[[pieces[1L]]] <- pieces[2L]
  }
  defaults
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1L])))

for (file in c(
  "00_config.R", "01_utils.R", "02_dgp.R", "03_estimators.R",
  "04_proposed.R", "05_metrics.R", "06_translasso.R", "07_scenarios.R",
  "08_reporting.R"
)) {
  source(file.path(script_dir, file), local = FALSE)
}

require_simulation_packages()
cli <- parse_cli(commandArgs(trailingOnly = TRUE))
simulation <- as.integer(cli$simulation)
cfg <- make_config(cli$mode)
seed_offset <- as.integer(cli$seed_offset)
if (is.na(seed_offset)) {
  stop("--seed_offset must be an integer.", call. = FALSE)
}
if (!cli$initial_tuning %in% c("joint", "separate")) {
  stop("--initial_tuning must be joint or separate.", call. = FALSE)
}
cfg$initial_tuning <- cli$initial_tuning
if (cli$grid == "fixed") {
  cfg$lambda_multipliers <- 1
  cfg$fixed_penalty_multipliers <- list(
    initial = 1,
    target = 1,
    pool = 1,
    lambda0 = 0.5,
    lambda1 = 0.5
  )
} else if (cli$grid == "narrow") {
  cfg$lambda_multipliers <- c(0.5, 1)
} else if (cli$grid == "compact") {
  cfg$initial_lambda_multipliers <- c(0.5, 1)
  cfg$target_lambda_multipliers <- c(0.5, 1)
  cfg$pool_lambda_multipliers <- c(0.5, 1, 2, 4, 8)
  cfg$transfer_lambda_multipliers <- c(0.5, 1, 2, 4, 8)
  cfg$linked_transfer_penalties <- TRUE
} else if (cli$grid == "wide") {
  cfg$lambda_multipliers <- c(0.25, 0.5, 1, 2, 4)
} else if (cli$grid == "extra_wide") {
  cfg$lambda_multipliers <- c(0.125, 0.25, 0.5, 1, 2, 4, 8)
} else if (cli$grid != "standard") {
  stop(paste(
    "--grid must be fixed, narrow, compact, standard, wide,",
    "or extra_wide."
  ),
       call. = FALSE)
}
if (!is.na(cli$replications)) {
  cfg$replications <- as.integer(cli$replications)
}
scenarios <- make_scenarios(simulation, cfg)
for (scenario_number in seq_along(scenarios)) {
  scenarios[[scenario_number]]$scenario_number <- scenario_number
}
if (cli$scenario != "all") {
  scenario_indices <- as.integer(
    strsplit(cli$scenario, ",", fixed = TRUE)[[1L]]
  )
  if (anyNA(scenario_indices) ||
      any(scenario_indices < 1L | scenario_indices > length(scenarios))) {
    stop("--scenario must be all or comma-separated valid scenario numbers.",
         call. = FALSE)
  }
  scenarios <- scenarios[scenario_indices]
}
methods <- methods_for_simulation(simulation)
fixed_parameters <- NULL
fixed_parameters_md5 <- NA_character_
if (!is.na(cli$fixed_params)) {
  if (!file.exists(cli$fixed_params)) {
    stop("Fixed-parameter file was not found: ", cli$fixed_params,
         call. = FALSE)
  }
  fixed_parameters <- data.table::fread(cli$fixed_params)
  required_fixed_columns <- c(
    "scenario_number", "parameter_role", "method", "source_index",
    "lambda_initial", "lambda0", "lambda1", "lambda_target",
    "trans_lambda_first", "trans_lambda_second"
  )
  missing_fixed_columns <- setdiff(
    required_fixed_columns, names(fixed_parameters)
  )
  if (length(missing_fixed_columns) > 0L) {
    stop(
      "Fixed-parameter file is missing columns: ",
      paste(missing_fixed_columns, collapse = ", "),
      call. = FALSE
    )
  }
  fixed_parameters_md5 <- unname(tools::md5sum(cli$fixed_params))
  cfg$implementation_version <- "2026-11-13-v8-pilot-fixed"
}
tuning_label <- if (is.null(fixed_parameters)) cli$grid else "pilot-fixed"
translasso_file <- if (is.na(cli$translasso)) NULL else cli$translasso
if ("TransLasso" %in% methods && is.null(translasso_file)) {
  stop(
    "Simulations 4 and 5 require --translasso=/path/to/TransLasso-functions.R",
    call. = FALSE
  )
}

workers <- as.integer(cli$workers)
if (is.na(workers) || workers < 1L) {
  stop("--workers must be a positive integer.", call. = FALSE)
}
workers <- min(workers, parallel::detectCores(logical = TRUE))
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

run_replication <- function(scenario, scenario_number, replication) {
  seed <- cfg$seed + seed_offset +
    scenario_number * 100000L + replication * 1000L
  local_cfg <- cfg
  local_cfg$K <- as.integer(scenario$K)
  local_cfg$scad_a <- 3.7 * local_cfg$K / 2
  local_cfg$n_fit <- as.integer(scenario$n_fit)
  local_cfg$n_source <- as.integer(scenario$n_source)
  if (!is.null(fixed_parameters)) {
    fixed_scenario_number <- as.integer(scenario_number)
    scenario_parameters <- fixed_parameters[
      fixed_parameters$scenario_number == fixed_scenario_number,
    ]
    method_parameters <- scenario_parameters[
      scenario_parameters$parameter_role == "method",
    ]
    source_parameters <- scenario_parameters[
      scenario_parameters$parameter_role == "source_initial",
    ]
    scenario_methods <- methods_for_scenario(simulation, scenario)
    required_methods <- setdiff(scenario_methods, "Proposed")
    if (!setequal(method_parameters$method, required_methods)) {
      stop(
        "Fixed method parameters do not match scenario ", scenario_number,
        call. = FALSE
      )
    }
    source_parameters <- source_parameters[
      order(source_parameters$source_index),
    ]
    if (nrow(source_parameters) != local_cfg$K ||
        !identical(as.integer(source_parameters$source_index),
                   seq_len(local_cfg$K))) {
      stop(
        "Fixed source-initial parameters do not match scenario ",
        scenario_number,
        call. = FALSE
      )
    }
    local_cfg$fixed_method_parameters <- method_parameters
    local_cfg$fixed_source_lambdas <- source_parameters$lambda_initial
  }
  local_cfg <- validate_config(local_cfg)

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
  validate_simulation_data(data, local_cfg)
  source_initials <- fit_source_initials(data$sources, local_cfg, seed + 10L)
  fitted <- list()
  scenario_methods <- methods_for_scenario(simulation, scenario)
  proposed_seed <- seed + 20L
  shared_tuning_seed <- proposed_seed + 300L

  if ("Proposed" %in% scenario_methods) {
    fitted$Proposed <- fit_proposed(
      data$target_fit, data$target_cal, data$sources, local_cfg,
      source_initials, proposed_seed
    )

    selected_method <- if (
      identical(fitted$Proposed$selected, "transfer")
    ) {
      "F-AdaTrans"
    } else {
      "Target-only"
    }
    if (selected_method %in% scenario_methods) {
      reused <- fitted$Proposed
      reused$model$method <- selected_method
      reused$selected <- NULL
      reused$r_transfer <- NULL
      reused$r_target <- NULL
      reused$outer_errors <- NULL
      fitted[[selected_method]] <- reused
    }
  }
  methods_to_fit <- setdiff(
    scenario_methods,
    c("Proposed", "TransLasso", names(fitted))
  )
  for (method in methods_to_fit) {
    fitted[[method]] <- fit_and_calibrate(
      method, data, local_cfg, source_initials, shared_tuning_seed
    )
  }
  if ("TransLasso" %in% scenario_methods) {
    fitted$TransLasso <- fit_and_calibrate_translasso(
      data, local_cfg, translasso_file, seed + 500L
    )
  }

  fitted <- fitted[scenario_methods[scenario_methods %in% names(fitted)]]

  result <- data.table::rbindlist(lapply(fitted, function(object) {
    evaluate_calibrated_method(object, data$target_test, data$truth$beta0)
  }))
  result[, `:=`(
    simulation = simulation,
    scenario_number = scenario_number,
    scenario = scenario$id,
    scenario_label = scenario$label,
    replication = replication,
    seed = seed,
    implementation_version = local_cfg$implementation_version,
    initial_tuning = local_cfg$initial_tuning,
    tuning_grid = tuning_label,
    fixed_parameters_md5 = fixed_parameters_md5
  )]
  result[]
}

run_suffix <- if (nzchar(cli$run_id)) paste0("_", cli$run_id) else ""
output_file <- file.path(
  script_dir,
  "results",
  sprintf("simulation_%s_%s%s.csv", simulation, cfg$mode, run_suffix)
)
if (file.exists(output_file) && tolower(cli$overwrite) != "true") {
  stop(
    "Output already exists: ", output_file,
    ". Use --overwrite=true to replace it.",
    call. = FALSE
  )
}

checkpoint_dir <- file.path(
  script_dir,
  "results",
  "checkpoints",
  sprintf("simulation_%s_%s%s", simulation, cfg$mode, run_suffix)
)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
resume <- tolower(cli$resume) == "true"

tasks <- do.call(rbind, lapply(seq_along(scenarios), function(scenario_number) {
  data.frame(
    scenario_index = scenario_number,
    scenario_number = scenarios[[scenario_number]]$scenario_number,
    replication = seq_len(cfg$replications)
  )
}))

run_task <- function(task_index) {
  scenario_index <- tasks$scenario_index[task_index]
  scenario_number <- tasks$scenario_number[task_index]
  replication <- tasks$replication[task_index]
  scenario <- scenarios[[scenario_index]]
  checkpoint_file <- file.path(
    checkpoint_dir,
    sprintf("%s_rep_%04d.rds", scenario$id, replication)
  )
  if (resume && file.exists(checkpoint_file)) {
    cached <- readRDS(checkpoint_file)
    expected_seed <- cfg$seed + seed_offset +
      scenario_number * 100000L + replication * 1000L
    compatible <- is.data.frame(cached) && nrow(cached) > 0L &&
      all(cached$simulation == simulation) &&
      all(cached$scenario == scenario$id) &&
      all(cached$replication == replication) &&
      all(cached$seed == expected_seed) &&
      "scenario_number" %in% names(cached) &&
      all(cached$scenario_number == scenario_number) &&
      "implementation_version" %in% names(cached) &&
      all(cached$implementation_version == cfg$implementation_version) &&
      "initial_tuning" %in% names(cached) &&
      all(cached$initial_tuning == cfg$initial_tuning) &&
      "tuning_grid" %in% names(cached) &&
      all(cached$tuning_grid == tuning_label) &&
      "fixed_parameters_md5" %in% names(cached) &&
      all((is.na(cached$fixed_parameters_md5) &
           is.na(fixed_parameters_md5)) |
          cached$fixed_parameters_md5 == fixed_parameters_md5)
    if (compatible) {
      message("Reusing ", basename(checkpoint_file))
      return(cached)
    }
    message("Recomputing incompatible checkpoint ", basename(checkpoint_file))
  }

  message(
    "Running ", scenario$id,
    ", replication ", replication, "/", cfg$replications
  )
  result <- run_replication(scenario, scenario_number, replication)
  temporary_file <- paste0(checkpoint_file, ".tmp-", Sys.getpid())
  saveRDS(result, temporary_file)
  if (!file.rename(temporary_file, checkpoint_file)) {
    stop("Could not finalize checkpoint: ", checkpoint_file, call. = FALSE)
  }
  message("Completed ", scenario$id, ", replication ", replication)
  result
}

task_indices <- seq_len(nrow(tasks))
if (workers == 1L) {
  all_results <- lapply(task_indices, run_task)
} else {
  all_results <- parallel::mclapply(
    task_indices,
    run_task,
    mc.cores = workers,
    mc.preschedule = FALSE
  )
}
failed <- vapply(all_results, inherits, logical(1), what = "try-error")
if (any(failed)) {
  stop(sum(failed), " parallel task(s) failed. Check the preceding log.",
       call. = FALSE)
}
results <- data.table::rbindlist(all_results, use.names = TRUE, fill = TRUE)

data.table::fwrite(results, output_file)
report <- write_simulation_outputs(
  results,
  output_stem = tools::file_path_sans_ext(output_file),
  alpha = cfg$alpha
)
print(results)
message("Saved results to ", output_file)
message("Saved tables to ", dirname(report$table_long))
message("Saved figures to ", report$figure_dir)
