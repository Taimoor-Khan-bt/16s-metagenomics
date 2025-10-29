# 16s-metagenomics – Release Notes

## v2.0.0 (2025-10-30)

Second public release: modular, config-driven pipeline with unified visualizations and cohort-scoped outputs.

Highlights
- Modular runner: `scripts/runner.R` orchestrates setup → preprocessing → analysis → visualization via `config/config.yaml`.
- Unified figures: single `scripts/visualization.R` controls all plots with `plots.enable` toggles (alpha, composition, heatmap, ordination, phylogenetic trees, per-sample panels).
- Cohort-scoped outputs: results saved under `output/<cohort>/{trimmed,filtered,analysis,visualizations}`; cohort defaults to `basename(io.input_dir)`.
- Optional tree building: DECIPHER alignment + phangorn NJ→ML (GTR+G+I) with circular/rectangular plots (ggtree).
- Warning hygiene: suppress non-actionable warnings; retain only upstream package deprecations.
- Documentation: README rewritten with quick start, configuration essentials, outputs, and biological interpretation.

Breaking changes
- Legacy monolithic scripts removed: `16s-metagenomics-complete.R`, `16s-metagenomics.R`, `16s_pipeline.yaml`.
- Visualization split replaced by unified `scripts/visualization.R`; old helpers removed.
- Output layout changed to cohort-scoped directories; scripts and config paths updated.

Upgrade notes
- Copy and edit `config/config.yaml` for your study; set `io.input_dir` and `amplicon.taxonomy` file paths.
- Use `plots.enable` toggles to customize figures; set `amplicon.phylogeny.build_tree: true` to enable tree building.
- Run via `Rscript scripts/runner.R --config config/config.yaml`. Optionally set `project.first_run: true` on new machines.

## v1.0.0 (2025-10-28)

Initial public release of the 16S rRNA metagenomics analysis pipeline (monolithic script).

Highlights
- Single entrypoint: `16s-metagenomics-complete.R`
- Robust cutadapt-based primer trimming
- DADA2 denoising with relaxed, cutadapt-aware filtering (no hard truncLen; tuned maxEE; minLen=100)
- SILVA taxonomy assignment (+ optional species)
- Phyloseq objects saved; alpha/beta diversity with statistical testing
- Publication-ready figures (TIFF, 600 dpi); rarefaction and length summaries
