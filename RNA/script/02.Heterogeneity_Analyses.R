# Figure 2c-2d
# =============================================================================
# Transcriptional similarity analysis: this study (Eu vs Aneu) vs Yan 2013 (Eu)
#   1) Load + QC + log2-transform each dataset
#   2) Take intersection of genes, merge, fix cell order
#   3) ComBat batch correction
#   4) HVG-based Pearson correlation
#   5) Heatmap + violin plots of correlation across stages
# =============================================================================

rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(sva)
  library(Matrix)
  library(matrixStats)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggpubr)
  library(pheatmap)
  library(RColorBrewer)
  library(grid)
})

# ---------- Config -----------------------------------------------------------
PATHS <- list(
  # datasets from this study: Smart-seq2 (81 cells) + STRT-seq (229 cells)，两份都是 raw counts
  this_counts_full    = "data/counts.matrix.full.QC.no_WYZX_E4.81.txt",
  this_counts_pooling = "data/counts.matrix.pooling.QC.with_WYZX_E4.229.txt",
  this_meta           = "data/metadata.310.txt",
  # Public dataset: Yan 2013 NSMB
  yan_expr            = "/Users/zhangjiaqi/Desktop/multi-omics/online_data/2013_NSMB_Yan/info/2013NSMB_count.txt",
  yan_meta            = "/Users/zhangjiaqi/Desktop/multi-omics/online_data/2013_NSMB_Yan/info/NSMB_sample_info.xls"
)

QC <- list(
  min_cells_per_gene = 3,      
  min_features       = 4000,   
  min_lib_size       = 1e5,    
  max_mito_frac      = 0.30,   
  mito_pattern       = "^MT-|^mt-"
)

STAGE_ORDER <- c("zygote", "2cell", "4cell", "8cell", "morula", "blastocyst")


YAN_KEEP_STAGES   <- c("h_4C", "h_8C")
YAN_EXCLUDE_CELLS <- c("h_8C_E2_6", "h_Mor_E1_3", "h_Mor_E1_8",
                       "h_Mor_E1_1", "h_Mor_E1_4") # Potentail aneuploid cells inferred by inferCNV, cut-off=0.8-1.2

# ---------- Helpers ----------------------------------------------------------
load_matrix <- function(path, sep = "\t") {
  df <- as.data.frame(fread(path, sep = sep, header = TRUE))
  rownames(df) <- df[[1]]
  df[[1]] <- NULL
  df
}

cell_qc <- function(counts,
                    min_features = 4000,
                    min_lib_size = 1e5,
                    max_mito     = 0.30,
                    mito_pattern = "^MT-|^mt-",
                    verbose      = TRUE) {
  counts <- as.matrix(counts)

  n_features <- Matrix::colSums(counts >= 1)
  lib_size   <- Matrix::colSums(counts)

  mito_genes <- grep(mito_pattern, rownames(counts), value = TRUE)
  if (length(mito_genes) == 0) {
    warning("cell_qc: no mitochondrial genes matched; skipping mito filter")
    mito_frac <- rep(0, ncol(counts))
  } else {
    mito_frac <- Matrix::colSums(counts[mito_genes, , drop = FALSE]) /
                 pmax(lib_size, 1)
  }

  pass_features <- n_features >= min_features
  pass_libsize  <- lib_size   >= min_lib_size
  pass_mito     <- mito_frac  <= max_mito
  keep <- pass_features & pass_libsize & pass_mito

  if (verbose) {
    cat(sprintf(
      "cell_qc: total=%d | fail n_features=%d | fail lib_size=%d | fail mito=%d | keep=%d\n",
      ncol(counts),
      sum(!pass_features),
      sum(!pass_libsize),
      sum(!pass_mito),
      sum(keep)
    ))
  }
  counts[, keep, drop = FALSE]
}

gene_filter_only <- function(mat, min_cells = 3) {
  keep <- Matrix::rowSums(mat >= 1) >= min_cells
  mat[keep, , drop = FALSE]
}

# preprocess pipeline: raw counts -> gene_filter -> CPM -> log2(CPM+1)
preprocess <- function(counts) {
  counts <- cell_qc(counts)
  counts <- gene_filter_only(counts)
  cpm_mat <- edgeR::cpm(counts, log = FALSE)
  log2(cpm_mat + 1)
}

# =============================================================================
# 1. This study --------------------------------------------------------------
counts_full    <- load_matrix(PATHS$this_counts_full)     # 81 cells
counts_pooling <- load_matrix(PATHS$this_counts_pooling)  # 229 cells
common_genes <- intersect(rownames(counts_full), rownames(counts_pooling))
data <- cbind(counts_full[common_genes, , drop = FALSE],
              counts_pooling[common_genes, , drop = FALSE])
data <- as.matrix(data)
cat("This study raw:", dim(data)[1], "genes x", dim(data)[2], "cells\n")
data <- preprocess(data)

meta <- load_matrix(PATHS$this_meta)
if ("ploidy" %in% colnames(meta) && !"state" %in% colnames(meta)) {
  meta$state <- meta$ploidy
}
meta <- meta[, c("stage", "embryo", "state"), drop = FALSE]
meta$batch <- "this study"

common_cells <- intersect(colnames(data), rownames(meta))
data <- data[, common_cells, drop = FALSE]
meta <- meta[common_cells, , drop = FALSE]

cat("This study after QC:", dim(data)[1], "genes x", dim(data)[2], "cells\n")

# =============================================================================
# 2. Public dataset (Yan 2013) -----------------------------------------------
# =============================================================================
yan_raw <- load_matrix(PATHS$yan_expr)
keep_cols <- grep(paste(YAN_KEEP_STAGES, collapse = "|"), colnames(yan_raw), value = TRUE)
keep_cols <- setdiff(keep_cols, YAN_EXCLUDE_CELLS)
yan_raw   <- yan_raw[, keep_cols, drop = FALSE]

online_data <- preprocess(yan_raw)

online_meta <- read.table(PATHS$yan_meta, row.names = 2, header = FALSE, sep = "\t")
online_meta <- online_meta[colnames(online_data), , drop = FALSE]
online_meta <- transform(
  online_meta,
  stage  = gsub("late_blastocyst", "blastocyst", V6),
  embryo = sub("_[^_]*$", "", rownames(online_meta)),
  state  = "Eu",
  batch  = "2013_NSMB"
)
online_meta <- online_meta[, c("stage", "embryo", "state", "batch"), drop = FALSE]

# =============================================================================
# 3. Merge ------------------------------
# =============================================================================
common_genes <- intersect(rownames(data), rownames(online_data))
combined <- cbind(data[common_genes, ], online_data[common_genes, ])

meta$sample        <- rownames(meta)
online_meta$sample <- rownames(online_meta)
meta_all <- rbind(meta, online_meta)
rownames(meta_all) <- meta_all$sample

meta_all$stage <- factor(meta_all$stage, levels = STAGE_ORDER)
meta_all       <- meta_all[order(meta_all$stage), , drop = FALSE]

combined <- combined[, rownames(meta_all), drop = FALSE]
stopifnot(identical(colnames(combined), rownames(meta_all)))
cat("Merged matrix:", dim(combined)[1], "genes x", dim(combined)[2], "cells\n")

# =============================================================================
# 4. ComBat batch correction --------------------------------------------------
# =============================================================================
mod          <- model.matrix(~ stage, data = meta_all)
combat_expr  <- ComBat(dat       = as.matrix(combined),
                       batch     = meta_all$batch,
                       mod       = mod,
                       par.prior = TRUE,
                       prior.plots = FALSE)

# =============================================================================
# 5. Pearson correlation on HVGs ---------------------------------------------
# =============================================================================
gene_sd  <- matrixStats::rowSds(combat_expr)
hvg      <- names(sort(gene_sd, decreasing = TRUE))[1:2000]
cor.val  <- cor(combat_expr[hvg, ], use = "pairwise.complete.obs", method = "pearson")

# =============================================================================
# 6. Heatmap ------------------------------------------------------------------
# =============================================================================
annotation_col <- data.frame(
  Stage  = factor(meta_all$stage, levels = STAGE_ORDER),
  Ploidy = factor(meta_all$state),
  Batch  = factor(meta_all$batch),
  row.names = rownames(meta_all)
)

ann_colors <- list(
  Stage  = c(zygote = "#E64B35", `2cell` = "#4DBBD5", `4cell` = "#00A087",
             `8cell` = "#3C5488", morula = "#F39B7F", blastocyst = "#8491B4"),
  Ploidy = c(Eu = "#1f77b4", Aneu = "#ff7f0e"),
  Batch  = c(`this study` = "#607D8B", `2013_NSMB` = "#CFD8DC")
)

pheatmap(
  cor.val,
  annotation_col    = annotation_col,
  annotation_colors = ann_colors,
  color             = colorRampPalette(rev(brewer.pal(7, "RdBu")))(100),
  breaks            = seq(-1, 1, length.out = 101),
  border_color      = NA,
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  fontsize          = 10,
  main              = "Pearson Correlation"
)

# =============================================================================
# 7. Transcriptional similarity: Eu vs Aneu against previous-stage Eu --------
# =============================================================================
group <- data.frame(
  cell   = colnames(combat_expr),
  stage  = meta_all$stage,
  embryo = meta_all$embryo,
  ploidy = meta_all$state,
  batch  = meta_all$batch,
  stringsAsFactors = FALSE
)

get_cells <- function(stage_val, ploidy_val) {
  group$cell[group$stage == stage_val & group$ploidy == ploidy_val]
}

build_pair_df <- function(query_cells, ref_cells, label, stage_val) {
  if (length(query_cells) == 0 || length(ref_cells) == 0) return(NULL)
  m  <- cor.val[query_cells, ref_cells, drop = FALSE]
  dt <- as.data.table(m, keep.rownames = "query_cell")
  dt <- melt(dt, id.vars = "query_cell",
             variable.name = "ref_cell", value.name = "correlation")
  dt$group <- label
  dt$stage <- stage_val
  dt
}

cor_df_list <- list()
for (i in 2:length(STAGE_ORDER)) {
  cur  <- STAGE_ORDER[i]
  prev <- STAGE_ORDER[i - 1]
  eu_prev <- get_cells(prev, "Eu")

  cor_df_list[[length(cor_df_list) + 1]] <-
    build_pair_df(get_cells(cur, "Eu"),   eu_prev, "Eu vs prev_Eu",   cur)
  cor_df_list[[length(cor_df_list) + 1]] <-
    build_pair_df(get_cells(cur, "Aneu"), eu_prev, "Aneu vs prev_Eu", cur)
}
cor_all <- rbindlist(Filter(Negate(is.null), cor_df_list))
cor_all[, group := factor(group, levels = c("Eu vs prev_Eu", "Aneu vs prev_Eu"))]
cor_all[, stage := factor(stage, levels = STAGE_ORDER)]

cor_plot <- cor_all[stage %in% c("8cell", "morula", "blastocyst")]

ggplot(cor_plot, aes(x = group, y = correlation, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.8, width = 0.9) +
  geom_boxplot(width = 0.2, outlier.shape = NA, color = "black",
               fill = NA, size = 0.4) +
  stat_summary(fun = median, geom = "point",
               shape = 95, size = 6, color = "black") +
  stat_compare_means(
    aes(label = after_stat(p.signif)),
    method      = "wilcox.test",
    label.y.npc = "top",
    label.x.npc = "center",
    size        = 3.5
  ) +
  facet_wrap(~ stage, nrow = 1, scales = "free_x") +
  scale_fill_manual(
    values = c("Eu vs prev_Eu" = "#1f77b4", "Aneu vs prev_Eu" = "#ff7f0e"),
    name   = "Comparison",
    labels = c("Eu vs previous-stage Eu", "Aneu vs previous-stage Eu")
  ) +
  labs(x = NULL, y = "Pearson correlation") +
  theme_classic(base_size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(size = 12),
    axis.text.x      = element_blank(),
    axis.text.y      = element_text(size = 12),
    axis.ticks.x     = element_blank(),
    legend.position  = "top",
    legend.title     = element_text(size = 9),
    legend.text      = element_text(size = 9),
    plot.margin      = margin(5, 10, 5, 5)
  )
