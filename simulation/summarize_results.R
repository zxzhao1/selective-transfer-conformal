script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[1L])))
source(file.path(script_dir, "08_reporting.R"), local = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop("Usage: Rscript summarize_results.R /path/to/results.csv [alpha]",
       call. = FALSE)
}
input_file <- normalizePath(args[1L])
alpha <- if (length(args) >= 2L) as.numeric(args[2L]) else 0.1
results <- data.table::fread(input_file)
report <- write_simulation_outputs(
  results,
  output_stem = tools::file_path_sans_ext(input_file),
  alpha = alpha
)
print(report$summary)
message("Saved tables to ", dirname(report$table_long))
message("Saved figures to ", report$figure_dir)
