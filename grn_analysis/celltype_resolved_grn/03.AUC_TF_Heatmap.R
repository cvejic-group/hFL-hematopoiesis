###########################################################
### Plot AUCell Heatmap with TF Expression + AUC scores ###
###########################################################

###-------------------Note-------------------###
###-------      Run with R-4.4.2      -------###
###------------------------------------------###

##################
# Initialization #
##################

setwd("~/work/")
source("./00.Initialization.R")

# Load data
TOPeReg_TF_AUCell_mat <- readRDS(paste0(RES_DIR, "TOPeReg_TF_AUCell_mat.rds"))
TF_EXPR.m <- TOPeReg_TF_AUCell_mat$TF_EXPR.m
AUC_GENE.m <- TOPeReg_TF_AUCell_mat$AUC_GENE.m
AUC_REGION.m <- TOPeReg_TF_AUCell_mat$AUC_REGION.m

# snRNA-seq data
FL.SeuratObj <- readRDS("~/work/FL_scrna_seurat.rds")
## Cell Metadata
cell_metadata <- FL.SeuratObj@meta.data

####################
# Data preparation #
####################

# Compute AVG
## EXPR
avg_expr <- TF_EXPR.m %>%
  t() %>%
  as.data.frame() %>%
  mutate(celltype = cell_metadata$anno_wnn_v51) %>%
  group_by(celltype) %>%
  summarise(across(everything(), mean)) %>%
  column_to_rownames(var = "celltype") %>%
  t()
zscore_expr <- t(apply(avg_expr, 1, function(x) {
  (x - mean(x)) / sd(x)
}))
norm_expr_df_TF <- as.data.frame(zscore_expr)
norm_expr_df_TF <- norm_expr_df_TF[, CT_ORDER]
## AUC_GENE
avg_expr <- AUC_GENE.m %>%
  t() %>%
  as.data.frame() %>%
  mutate(celltype = cell_metadata$anno_wnn_v51) %>%
  group_by(celltype) %>%
  summarise(across(everything(), mean)) %>%
  column_to_rownames(var = "celltype") %>%
  t()  
### Min-Max norm
zscore_expr <- t(apply(avg_expr, 1, function(x) {
  (x - mean(x)) / sd(x)
}))
norm_expr <- t(apply(zscore_expr, 1, function(x) {
  (x - min(x)) / (max(x) - min(x))
}))
norm_expr_df_gene <- as.data.frame(norm_expr)
norm_expr_df_gene <- norm_expr_df_gene[, CT_ORDER]
## AUC_REGION
avg_expr <- AUC_REGION.m %>%
  t() %>%
  as.data.frame() %>%
  mutate(celltype = cell_metadata$anno_wnn_v51) %>%
  group_by(celltype) %>%
  summarise(across(everything(), mean)) %>%
  column_to_rownames(var = "celltype") %>%
  t()  
### Min-Max norm
zscore_expr <- t(apply(avg_expr, 1, function(x) {
  (x - mean(x)) / sd(x)
}))
norm_expr <- t(apply(zscore_expr, 1, function(x) {
  (x - min(x)) / (max(x) - min(x))
}))
norm_expr_df_region <- as.data.frame(norm_expr)
norm_expr_df_region <- norm_expr_df_region[, CT_ORDER]

# Get Gene Order
norm_expr_df_gene2 <- as.matrix(norm_expr_df_gene)
norm_expr_df_gene2[which(norm_expr_df_gene2 < 1)] <- 0
gene_order <- c()
for (i in 1:ncol(norm_expr_df_gene2)) {
  gene_order <- 
    names(which(norm_expr_df_gene2[, i] > 0)) %>%
    c(gene_order, .) 
}
tmp <- as.matrix(norm_expr_df_TF)
## Set EXPR Z-score boundary
tmp[which(tmp < -2)] <- -2
tmp[which(tmp > 2)] <- 2

###################
# Plot Generation #
###################

# Prepare plot data
df_color <- as.data.frame(tmp, check.names = FALSE) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "CellType", values_to = "TF")
df_size <- as.data.frame(norm_expr_df_gene, check.names = FALSE) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "CellType", values_to = "Gene_base")
df_alpha <- as.data.frame(norm_expr_df_region, check.names = FALSE) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "CellType", values_to = "Region_base")
df_all <- purrr::reduce(list(df_color, df_size, df_alpha), dplyr::left_join, by = c("Gene", "CellType"))

# Generate Plots
## With TF Expression
Cairo::CairoPDF(paste0(FIG_DIR, "FL_AUC_TF_heatmap.pdf"), 
                width = 48, height = 12, family = "Arial")
ggplot(df_all, aes(x = CellType, y = Gene)) +
  geom_tile(aes(fill = TF), color = "white") +
  geom_point(aes(size = Gene_base, alpha = Region_base), color = "black") +
  scale_fill_gradientn(
    colours = c("#3a86ff", "#48cae4", "white", "#e9c46a", "#d62828"),
    values  = scales::rescale(c(-2, -1, 0, 1, 2)),
    name    = "TF expression"
  ) +
  scale_size(range = c(1, 5), name = "Gene-based AUCell scores") +
  scale_alpha(range = c(0.1, 1), name = "Region-based AUCell scores") +
  scale_y_discrete(limits = gene_order) +
  scale_x_discrete(limits = rev(CT_ORDER)) +
  labs(y = "TF", x = "Cell Types", title = "AUCell Scores across Cell Types") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = .5, colour = "black", size = 17)
  ) +
  coord_flip()
dev.off()
## Without TF Expression
Cairo::CairoPDF(paste0(FIG_DIR, "FL_AUC_heatmap.pdf"), 
                width = 48, height = 12, family = "Arial")
ggplot(df_all, aes(x = as.factor(CellType), y = Gene)) +
  geom_tile(aes(fill = Gene_base), color = "white") +
  geom_point(aes(size = Region_base), color = "black") +
  scale_fill_gradientn(
    colours = c("#3a86ff", "#48cae4", "white", "#e9c46a", "#d62828"),
    values = scales::rescale(c(0, .25, .5, .75, 1)),
    name = "Gene-based AUCell scores"
  ) +
  scale_size(range = c(1, 4), name = "Region-based AUCell scores") +
  scale_y_discrete(limits = gene_order) +
  scale_x_discrete(limits = rev(CT_ORDER)) +
  theme_minimal(base_size = 12) +
  labs(y = "TF", x = "Cell Types", title = "AUCell Scores across Cell Types") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(face = "bold"),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = .5, colour = "black", size = 17)) +
  coord_flip()
dev.off()


