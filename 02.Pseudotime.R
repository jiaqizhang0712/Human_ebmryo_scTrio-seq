# 3. Monocle3 Object
rm(list = ls())

library(monocle3)
packageVersion("monocle3") # 1.3.7
library(Seurat)
packageVersion("Seurat") # 5.3.0
library(dplyr)
library(ggplot2)
library(patchwork)
library(Matrix)

counts1 <- read.table("data/counts.matrix.full.QC.no_WYZX_E4.81.txt", row.names = 1, header = TRUE, sep = "\t") 
counts2 <- read.table("data/counts.matrix.pooling.QC.with_WYZX_E4.229.txt", row.names = 1, header = TRUE, sep = "\t")
dim(counts1) # 43901   81
dim(counts2) # 59505   229
common_genes <- intersect(rownames(counts1), rownames(counts2))
length(common_genes) # 43901
exprs <- cbind(counts1[common_genes, , drop = FALSE], counts2[common_genes,  , drop = FALSE])
dim(exprs) #43901   310
exprs <- exprs[rowSums(exprs > 0) > 0, ]
dim(exprs) # 40401   310

pheno_data <- read.table("data/metadata.310.txt", row.names = 1, header = TRUE, sep = "\t", check.names = FALSE)
pheno_data <- pheno_data[rownames(pheno_data) %in% colnames(exprs), ]
exprs <- exprs[, rownames(pheno_data), drop = FALSE]

gene_annotation <- data.frame(gene_short_name = rownames(exprs),row.names = rownames(exprs))

cds <- new_cell_data_set(
  as(as.matrix(exprs), "sparseMatrix"),
  cell_metadata = pheno_data,
  gene_metadata = gene_annotation
) 

# 4. Trajectory 
set.seed(123) 
cds <- preprocess_cds(cds, num_dim = 15)
cds <- align_cds(cds, num_dim = 10, alignment_group = "library") # batch "mutual nearest neighbor" algorithm PMID: 29608177
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds, reduction_method = "UMAP", resolution = 1:10) 
cds <- learn_graph(cds)
cds <- order_cells(cds, root_cells = "ZHYZ_E1") # zygote

# 5. Figure 2a
p1 <- plot_cells(cds, 
                 color_cells_by="stage",
                 group_label_size = 2.5,
                 cell_size = 1,
                 label_cell_groups = TRUE,
                 label_groups_by_cluster=TRUE,
                 label_leaves=TRUE,
                 label_branch_points=TRUE,
                 graph_label_size=1.5,
                 show_trajectory_graph = TRUE,)

p2 <- plot_cells(cds, 
                 color_cells_by="pseudotime",
                 group_label_size = 2.5,cell_size = 1,
                 label_cell_groups = TRUE,
                 label_groups_by_cluster=TRUE,
                 label_leaves=TRUE,
                 label_branch_points=TRUE,
                 graph_label_size=1.5,
                 show_trajectory_graph = TRUE,)

print(p1|p2)


# results
pseudotime <- pseudotime(cds, reduction_method = "UMAP")
cds@colData$Pseudotime <- pseudotime
cell_meta <- as.data.frame(cds@colData)
cell_meta$cell <- rownames(cell_meta)
cell_meta_f <- cell_meta
plot_df <- cell_meta_f[, c("cell", "stage", "ploidy", "Pseudotime")]
plot_df$stage  <- factor(plot_df$stage, levels = c("zygote","2cell","4cell","8cell","morula","blastocyst"))
plot_df$ploidy <- factor(plot_df$ploidy, levels = c("Eu","Aneu"))
df_zygote <- dplyr::filter(plot_df, stage == "zygote")
df_2cell  <- dplyr::filter(plot_df, stage == "2cell")
df_4cell  <- dplyr::filter(plot_df, stage == "4cell")
df_8cell  <- dplyr::filter(plot_df, stage == "8cell")
df_morula <- dplyr::filter(plot_df, stage == "morula")
df_blast  <- dplyr::filter(plot_df, stage == "blastocyst")

# Figure 2b
library(ggunchained) 
if(T){
  mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                   axis.title = element_text(size = 12,color ="black"), 
                   axis.text = element_text(size= 12,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(angle = 0, hjust = 1 ),
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12)
  ) 
}
ggplot(plot_df, aes(x = stage, y = Pseudotime, fill = ploidy)) +
  
  geom_split_violin(data = df_zygote, trim = FALSE, color = "grey", scale = "area") +
  geom_split_violin(data = df_2cell,  trim = FALSE, color = "grey", scale = "area") +
  geom_split_violin(data = df_4cell,  trim = FALSE, color = "grey", scale = "area") +
  geom_split_violin(data = df_8cell,  trim = FALSE, color = "grey", scale = "area") +
  geom_split_violin(data = df_morula, trim = FALSE, color = "grey", scale = "area") +
  geom_split_violin(data = df_blast,  trim = FALSE, color = "grey", scale = "area") +
  
  # median points
  geom_point(
    data = plot_df,
    aes(x = stage, y = Pseudotime, fill = ploidy),
    stat = "summary",
    fun = median,
    position = position_dodge(width = 0.2)
  ) +
  
  # IQR error bars (Q1–Q3)
  stat_summary(
    data = plot_df,
    aes(x = stage, y = Pseudotime, group = ploidy),
    fun.min = function(x) quantile(x)[2],
    fun.max = function(x) quantile(x)[4],
    geom = "errorbar",
    color = "black",
    width = 0.1,
    size = 0.5,
    position = position_dodge(width = 0.2)
  ) +
  
  # Wilcoxon test per stage (Eu vs Aneu)
  stat_compare_means(
    aes(x = stage, y = Pseudotime, group = ploidy),
    method = "wilcox",
    label = "p.signif",
    label.y = max(plot_df$Pseudotime, na.rm = TRUE) + 5,
    hide.ns = TRUE,
    symnum.args = list(
      cutpoints = c(0, 0.001, 0.01, 0.05, 1),
      symbols = c("***", "**", "*", "NS")
    ),
    size = 3.5
  ) +
  
  scale_fill_manual(values = c("#56B4E9", "#E69F00")) +
  theme_bw() +
  labs(y = "Pseudotime", x = "Stage") +
  mytheme +
  theme(
    legend.position = c(0.1, 0.02),
    legend.justification = c("left", "bottom"),
    panel.border = element_rect(color = "black", size = 0.8)
  ) +
  guides(fill = guide_legend(override.aes = list(shape = 20, size = 2))) +
  coord_flip()
