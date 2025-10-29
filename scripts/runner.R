# Orchestrator for metagenomics analyses (16S | shotgun)
# Usage: Rscript scripts/runner.R --config config/config.yaml

suppressPackageStartupMessages({
  library(yaml)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || args[1] != "--config") {
  stop("Usage: Rscript scripts/runner.R --config <config.yaml>")
}
cfg_path <- normalizePath(args[2])
cfg <- yaml::read_yaml(cfg_path)

message("[runner] Project: ", cfg$project$name)
message("[runner] Mode: ", cfg$project$sequencing_type)

dir.create(cfg$project$output_dir, showWarnings = FALSE, recursive = TRUE)
if (!is.null(cfg$project$temp_dir)) dir.create(cfg$project$temp_dir, showWarnings = FALSE, recursive = TRUE)

`%||%` <- function(a, b) if (!is.null(a)) a else b
set.seed(cfg$project$random_seed %||% 1234)

if (tolower(cfg$project$sequencing_type) == "16s") {
  # Setup if requested
  if (isTRUE(cfg$project$first_run)) {
    source(file.path("scripts", "setup.R"))
    run_setup(cfg)
  } else {
    # Ensure basic cohort directories exist
    cohort <- cfg$io$cohort; if (is.null(cohort) || cohort == "") cohort <- basename(cfg$io$input_dir)
    out_base <- file.path(cfg$project$output_dir, cohort)
    dir.create(file.path(out_base, "trimmed"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(out_base, "filtered"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(out_base, "analysis"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(out_base, "visualizations"), recursive = TRUE, showWarnings = FALSE)
  }

  # Pipeline stages
  source(file.path("scripts", "preprocess_16s.R"))
  source(file.path("scripts", "analysis_16s.R"))
  source(file.path("scripts", "visualization.R"))

  pre <- run_preprocess_16s(cfg)
  run_analysis_16s(cfg, pre)
  run_plots(cfg)

} else if (tolower(cfg$project$sequencing_type) == "shotgun") {
  # Execute shotgun module (bash) then visualize in R
  sh <- file.path("scripts", "module_shotgun.sh")
  if (!file.exists(sh)) stop("Missing scripts/module_shotgun.sh")

  # Pass key config as CLI args to the bash module
  bash_args <- c(
    sh,
    "--input-dir", cfg$io$input_dir,
    "--output-dir", cfg$project$output_dir,
    "--threads", as.character(cfg$project$threads %||% 4),
    "--metadata", cfg$io$metadata_csv,
    "--host-ref", cfg$shotgun$host_depletion$reference %||% "GRCh38",
    "--taxonomy", cfg$shotgun$taxonomy$tool %||% "metaphlan4",
  "--func-tool", cfg$shotgun[["function"]][["tool"]] %||% "humann",
    "--min-read-len", as.character(cfg$shotgun$min_read_len %||% 50)
  )
  message("[runner] Launching shotgun module: ", paste(bash_args, collapse = " "))
  status <- system2("bash", bash_args, stdout = TRUE, stderr = TRUE)
  cat(paste(status, collapse = "\n"), "\n")

  source(file.path("scripts", "visualization_common.R"))
  source(file.path("scripts", "visualization_shotgun.R"))
  run_shotgun_plots(cfg)

} else {
  stop("Unsupported sequencing_type in config: ", cfg$project$sequencing_type)
}

message("[runner] Done.")
