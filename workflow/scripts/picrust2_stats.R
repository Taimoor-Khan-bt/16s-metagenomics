#!/usr/bin/env Rscript
# =============================================================================
# picrust2_stats.R — Robust Multi-group Differential Pathway Analysis
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# ── Arguments ─────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: picrust2_stats.R <pathways_tsv> <metadata> <group_col> <out_dir>")

pathways_file <- args[1]
meta_file     <- args[2]
group_col     <- args[3]
out_dir       <- args[4]

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── 1. Load Data ─────────────────────────────────────────────────────────────
tab <- read.table(pathways_file, header = TRUE, sep = "\t",
                  row.names = 1, check.names = FALSE, comment.char = "#")

meta <- read.table(meta_file, header = TRUE, sep = "\t",
                   comment.char = "", row.names = 1, check.names = FALSE)

# Sync samples
common_samples <- intersect(colnames(tab), rownames(meta))
if (length(common_samples) < 3) stop("Fewer than 3 common samples. Check Sample IDs.")
tab  <- tab[, common_samples, drop = FALSE]
meta <- meta[common_samples, , drop = FALSE]

# Ensure group is a factor and remove levels with zero samples
meta[[group_col]] <- as.factor(meta[[group_col]])
grp  <- meta[[group_col]]
lvls <- levels(grp)

message("Analysis Groups: ", paste(lvls, collapse=", "))

# ── 2. Kruskal-Wallis Test (Omnibus) ─────────────────────────────────────────
# Correct for 3+ groups instead of just Wilcoxon
results <- lapply(rownames(tab), function(pw) {
  vals <- as.numeric(tab[pw, ])
  if (sum(vals) == 0) return(NULL)
  
  # Omnibus test: Is there any difference across the 3 groups?
  p_val <- tryCatch(kruskal.test(vals ~ grp)$p.value, error = function(e) NA)
  
  # Calculate means per group
  group_means <- tapply(vals, grp, mean)
  
  res <- data.frame(
    pathway = pw,
    p_value = p_val,
    stringsAsFactors = FALSE
  )
  # Append means dynamically for any number of groups
  for(l in lvls) { res[[paste0("mean_", l)]] <- group_means[l] }
  return(res)
})

res_df <- do.call(rbind, Filter(Negate(is.null), results))
res_df$p_adj <- p.adjust(res_df$p_value, method = "BH")
res_df <- res_df[order(res_df$p_adj), ]

write.table(res_df, file.path(out_dir, "pathway_differential.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ── 3. Visualization ──────────────────────────────────────────────────────────
# Select top 20 significant pathways
sig <- res_df %>% filter(!is.na(p_adj), p_adj < 0.25) %>% head(20)
if (nrow(sig) == 0) sig <- head(res_df, 20)

# Prepare data for plotting (Pivot to long format)
plot_df <- sig %>%
  select(pathway, starts_with("mean_")) %>%
  pivot_longer(cols = starts_with("mean_"), names_to = "Group", values_to = "Abundance") %>%
  mutate(Group = gsub("mean_", "", Group))

pdf(file.path(out_dir, "pathway_plots.pdf"), width = 14, height = 8)

# Grouped Bar Plot (Better for 3 groups than a Volcano plot)
p_bar <- ggplot(plot_df, aes(x = reorder(pathway, Abundance), y = Abundance, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  scale_fill_brewer(palette = "Set1") +
  theme_bw(base_size = 12) +
  labs(title = "Top Predicted Pathways (Abundance by Group)",
       subtitle = "Omnibus Kruskal-Wallis Test",
       x = "MetaCyc Pathway", y = "Relative Abundance") +
  theme(legend.position = "bottom")

print(p_bar)

# Heatmap of top pathways
p_heat <- ggplot(plot_df, aes(x = Group, y = pathway, fill = log10(Abundance + 1))) +
  geom_tile() +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(title = "Pathway Abundance Heatmap (Log10 scale)",
       x = "Group", y = "Pathway")

print(p_heat)

dev.off()
message("Saved: pathway_differential.tsv and pathway_plots.pdf")