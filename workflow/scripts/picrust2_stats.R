#!/usr/bin/env Rscript
# =============================================================================
# picrust2_stats.R — Differential pathway analysis on PICRUSt2 outputs
# =============================================================================
# Usage: Rscript picrust2_stats.R <path_abun_tsv> <metadata_tsv>
#                                 <group_col> <out_dir>
#
# Required: ggplot2, dplyr, tidyr
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

# ── Load data ─────────────────────────────────────────────────────────────────
tab <- read.table(pathways_file, header = TRUE, sep = "\t",
                  row.names = 1, check.names = FALSE, comment.char = "#")

meta <- read.table(meta_file, header = TRUE, sep = "\t",
                   comment.char = "", row.names = 1, check.names = FALSE)

# Tab: rows = pathways, cols = samples
common_samples <- intersect(colnames(tab), rownames(meta))
if (length(common_samples) < 2) stop("Fewer than 2 common samples.")
tab  <- tab[, common_samples, drop = FALSE]
meta <- meta[common_samples, , drop = FALSE]

grp    <- as.factor(meta[[group_col]])
lvls   <- levels(grp)

message("Samples: ", length(common_samples), " | Pathways: ", nrow(tab))

# ── Wilcoxon rank-sum per pathway ─────────────────────────────────────────────
g1_idx <- which(grp == lvls[1])
g2_idx <- which(grp == lvls[2])

results <- lapply(rownames(tab), function(pw) {
  x1 <- as.numeric(tab[pw, g1_idx])
  x2 <- as.numeric(tab[pw, g2_idx])
  if (all(x1 == 0) && all(x2 == 0)) return(NULL)
  p  <- tryCatch(wilcox.test(x1, x2)$p.value, error = function(e) NA)
  m1 <- mean(x1 + 1e-10); m2 <- mean(x2 + 1e-10)
  data.frame(
    pathway    = pw,
    log2fc     = log2(m1) - log2(m2),
    mean_group1 = mean(x1),
    mean_group2 = mean(x2),
    p_value    = p,
    stringsAsFactors = FALSE
  )
})

res_df <- do.call(rbind, Filter(Negate(is.null), results))
res_df$p_adj <- p.adjust(res_df$p_value, method = "BH")
res_df <- res_df[order(res_df$p_adj), ]
colnames(res_df)[3:4] <- paste0("mean_", lvls)

write.table(res_df, file.path(out_dir, "pathway_differential.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("Saved: pathway_differential.tsv (", nrow(res_df), " pathways)")

# ── Plots ─────────────────────────────────────────────────────────────────────
sig <- res_df %>% filter(!is.na(p_adj), p_adj < 0.25) %>% arrange(desc(abs(log2fc)))
if (nrow(sig) == 0) sig <- head(res_df, 30)

pdf(file.path(out_dir, "pathway_plots.pdf"), width = 12, height = max(6, nrow(sig) * 0.3 + 2))

# Volcano plot
p_volcano <- ggplot(res_df, aes(x = log2fc, y = -log10(p_value),
                                color = p_adj < 0.25)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("grey60", "tomato"),
                     labels = c("p_adj ≥ 0.25", "p_adj < 0.25")) +
  theme_bw(base_size = 12) +
  labs(title = paste("Differential Pathways —", lvls[1], "vs", lvls[2]),
       x = "Log2 Fold-Change", y = "-log10(p-value)", color = NULL) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey40") +
  geom_vline(xintercept = 0, linetype = 1, color = "grey40")
print(p_volcano)

# Top pathways barplot
p_bar <- ggplot(sig, aes(x = reorder(pathway, log2fc), y = log2fc,
                          fill = ifelse(log2fc > 0, lvls[1], lvls[2]))) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("steelblue", "tomato")) +
  theme_bw(base_size = 10) +
  theme(legend.title = element_blank()) +
  labs(title = paste("Top Differential Pathways (p_adj < 0.25)"),
       x = "Pathway", y = "Log2 Fold-Change")
print(p_bar)

dev.off()
message("Saved: pathway_plots.pdf")
