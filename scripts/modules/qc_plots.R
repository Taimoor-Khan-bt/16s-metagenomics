# scripts/modules/qc_plots.R

suppressPackageStartupMessages({
  library(fastqcr)
  library(ggplot2)
  library(patchwork)
})

generate_qc_plots <- function(cfg) {
  message("--- Generating Static QC Plots ---")

  output_base <- cfg$project$output_dir %||% "output"
  cohort <- cfg$io$cohort %||% basename(cfg$io$input_dir)
  outdir <- file.path(output_base, cohort, "visualizations", "qc")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  fastqc_raw_dir <- file.path(output_base, "fastqc_raw")
  fastqc_trimmed_dir <- file.path(output_base, "fastqc_trimmed")

  # Aggregate FastQC reports
  qc_raw <- qc_aggregate(fastqc_raw_dir)
  qc_trimmed <- qc_aggregate(fastqc_trimmed_dir)

  # --- Per Base Sequence Quality ---
  p1 <- qc_plot(qc_raw, "Per base sequence quality") +
    ggtitle("Raw Reads") +
    theme(legend.position = "none")
  p2 <- qc_plot(qc_trimmed, "Per base sequence quality") +
    ggtitle("Trimmed Reads") +
    theme(legend.position = "none")

  p_qual <- p1 + p2 + plot_layout(ncol = 2)
  ggsave(file.path(outdir, "per_base_quality.png"), p_qual, width = 12, height = 6, dpi = 300)

  # --- Per Sequence GC Content ---
  p1 <- qc_plot(qc_raw, "Per sequence GC content") +
    ggtitle("Raw Reads") +
    theme(legend.position = "none")
  p2 <- qc_plot(qc_trimmed, "Per sequence GC content") +
    ggtitle("Trimmed Reads") +
    theme(legend.position = "none")

  p_gc <- p1 + p2 + plot_layout(ncol = 2)
  ggsave(file.path(outdir, "per_sequence_gc.png"), p_gc, width = 12, height = 6, dpi = 300)

  # --- Sequence Length Distribution ---
  p1 <- qc_plot(qc_raw, "Sequence Length Distribution") +
    ggtitle("Raw Reads") +
    theme(legend.position = "none")
  p2 <- qc_plot(qc_trimmed, "Sequence Length Distribution") +
    ggtitle("Trimmed Reads") +
    theme(legend.position = "none")

  p_len <- p1 + p2 + plot_layout(ncol = 2)
  ggsave(file.path(outdir, "sequence_length_distribution.png"), p_len, width = 12, height = 6, dpi = 300)

  message("--- Static QC Plots Generated ---")
}
