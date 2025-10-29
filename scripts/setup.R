# Setup script: installs R packages, downloads reference files, and scaffolds directories

suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
  library(utils)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

install_if_missing <- function(pkgs_cran = character(), pkgs_bioc = character()) {
  quietly <- function(expr) suppressPackageStartupMessages(suppressWarnings(expr))
  for (p in pkgs_cran) {
    if (!requireNamespace(p, quietly = TRUE)) {
      message("[setup] Installing CRAN package: ", p)
      quietly(install.packages(p, repos = "https://cloud.r-project.org"))
    }
  }
  for (p in pkgs_bioc) {
    if (!requireNamespace(p, quietly = TRUE)) {
      message("[setup] Installing Bioconductor package: ", p)
      quietly(BiocManager::install(p, ask = FALSE, update = FALSE))
    }
  }
}

scaffold_dirs <- function(cfg) {
  out_base <- cfg$project$output_dir %||% "output"
  cohort <- cfg$io$cohort
  if (is.null(cohort) || is.na(cohort) || cohort == "") cohort <- basename(cfg$io$input_dir)
  out_cohort <- file.path(out_base, cohort)
  dirs <- c(out_cohort,
            file.path(out_cohort, "trimmed"),
            file.path(out_cohort, "filtered"),
            file.path(out_cohort, "analysis"),
            file.path(out_cohort, "visualizations"),
            file.path(out_cohort, "logs"))
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  if (!is.null(cfg$project$temp_dir)) dir.create(cfg$project$temp_dir, recursive = TRUE, showWarnings = FALSE)
  invisible(out_cohort)
}

download_if_missing <- function(url, dest) {
  if (is.null(url) || is.na(url) || url == "") return(invisible(FALSE))
  if (!file.exists(dest)) {
    message("[setup] Downloading ", basename(dest), " from ", url)
    try(utils::download.file(url, destfile = dest, mode = "wb", quiet = TRUE), silent = TRUE)
  }
}

run_setup <- function(cfg) {
  message("[setup] Starting first-run setup...")
  # 1) Install required packages
  if (isTRUE(cfg$setup$install_packages)) {
    install_if_missing(
      pkgs_cran = c("tidyverse","ggplot2","data.table","gridExtra","RColorBrewer","viridisLite","yaml","remotes","reshape2","scales"),
      pkgs_bioc = c("dada2","phyloseq","vegan","DECIPHER","phangorn","ggtree","ComplexHeatmap","indicspecies")
    )
  }
  # 2) Create directory scaffold
  out_cohort <- scaffold_dirs(cfg)
  # 3) Download taxonomy assets if URLs provided and files missing
  silva_train <- cfg$amplicon$taxonomy$silva_train_set %||% "silva_nr99_v138.1_train_set.fa"
  silva_species <- cfg$amplicon$taxonomy$silva_species %||% "silva_species_assignment_v138.1.fa"
  download_if_missing(cfg$setup$taxonomy_urls$silva_train_set %||% NULL, silva_train)
  download_if_missing(cfg$setup$taxonomy_urls$silva_species %||% NULL, silva_species)
  # 4) Optionally pull scripts from GitHub
  repo <- cfg$setup$scripts_repo %||% NULL
  if (!is.null(repo) && repo != "") {
    message("[setup] Pulling scripts from ", repo)
    # best-effort; skip if git missing or offline
    try(system2("git", c("clone", repo, "scripts_repo_tmp"), stdout = TRUE, stderr = TRUE), silent = TRUE)
  }
  message("[setup] Setup completed. Cohort output at: ", out_cohort)
}
