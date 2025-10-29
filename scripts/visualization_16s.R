## Deprecated: This file has been replaced by `scripts/visualization.R`.
## The unified visualization script consolidates helpers and plots with
## config-based enable/disable toggles. This stub remains only for backward
## compatibility in case external callers source it directly.

run_16s_plots <- function(cfg) {
  message("[viz] visualization_16s.R is deprecated. Delegating to scripts/visualization.R")
  source(file.path("scripts", "visualization.R"))
  run_plots(cfg)
}
