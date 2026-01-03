#!/bin/bash
set -e

# verify R is installed
if ! command -v Rscript &> /dev/null; then
    echo "Error: Rscript is not installed or not in PATH."
    exit 1
fi

echo "Starting R package installation..."
echo "This may take a while depending on the number of packages to compile."

Rscript -e '
# Function to install packages if missing
install_pkgs <- function(pkgs) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        message("Installing BiocManager...")
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    
    # Configure BiocManager to use binary repositories if possible (faster on Linux)
    # This is often handled by R configuration, but good to ensure valid repos are set.
    
    installed <- rownames(installed.packages())
    to_install <- pkgs[!pkgs %in% installed]
    
    if (length(to_install) > 0) {
        message("Installing missing packages: ", paste(to_install, collapse = ", "))
        BiocManager::install(to_install, update = FALSE, ask = FALSE)
    } else {
        message("All packages are already installed.")
    }
}

# List of packages used in the project
# Combined from scripts/setup.R and grep analysis of library() calls
required_packages <- c(
    # CRAN Packages
    "tidyverse",
    "data.table",
    "gridExtra",
    "RColorBrewer",
    "viridisLite",
    "yaml",
    "remotes",
    "reshape2",
    "scales",
    "fastqcr",
    "ggrepel",
    "optparse",
    "patchwork",
    "ape",
    "vegan",
    
    # Bioconductor Packages
    "dada2",
    "phyloseq",
    "DECIPHER",
    "phangorn",
    "ggtree",
    "ComplexHeatmap",
    "indicspecies",
    "Biostrings",
    "DESeq2"
)

install_pkgs(required_packages)
'

echo "Installation complete."
