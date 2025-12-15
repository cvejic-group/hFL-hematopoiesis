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

# Select for specific purposes
## Plot eRegulon activity for PCW5
### Cell selections
SUBSET = 5
CELL_NUM = 30
FL.subset.SeuratObj <- subset(FL.SeuratObj, 
                              subset = PCW == SUBSET)
CT.keep <- intersect(
  CT_ORDER,
  names(which(table(FL.subset.SeuratObj$celltype_latest) >= CELL_NUM))
)
##### Alternative: manually set CT
CT.keep <- c("HSC", "GP",
             "MEMP-t", "MastP-t",
             "MDP",
             "LMPP", "LP", "Cycling-LP", "PreProB")
FL.subset.SeuratObj <- subset(FL.subset.SeuratObj, 
                              subset = celltype_latest %in% CT.keep)
### TF selections
XGBoost_RMT_TFRes.l <- read_rds(paste0(RES_DIR, 
                                       "XGBoost_RMT_TFRes_wRankedFeature.rds"))
TF.v <- lapply(XGBoost_RMT_TFRes.l[CT.keep], 
               function(x){
                 return(x$ranked_features[1:5])}) |> 
  unlist() |> 
  unique()
### Coordinate
TF_EXPR_sub.m <- TF_EXPR.m[TF.v, colnames(FL.subset.SeuratObj)]
AUC_GENE_sub.m <- AUC_GENE.m[TF.v, colnames(FL.subset.SeuratObj)]
AUC_REGION_sub.m <- AUC_GENE.m[TF.v, colnames(FL.subset.SeuratObj)]
cell_metadata_sub <- FL.subset.SeuratObj@meta.data

# Compute AVG
## EXPR
avg_expr <- TF_EXPR_sub.m %>%
  t() %>%
  as.data.frame() %>%
  mutate(celltype = cell_metadata_sub$anno_wnn_v51) %>%
  group_by(celltype) %>%
  summarise(across(everything(), mean)) %>%
  column_to_rownames(var = "celltype") %>%
  t()
zscore_expr <- t(apply(avg_expr, 1, function(x) {
  (x - mean(x)) / sd(x)
}))
norm_expr_df_TF <- as.data.frame(zscore_expr)
norm_expr_df_TF <- norm_expr_df_TF[, CT.keep]
## AUC_GENE
avg_expr <- AUC_GENE_sub.m %>%
  t() %>%
  as.data.frame() %>%
  mutate(celltype = cell_metadata_sub$anno_wnn_v51) %>%
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
norm_expr_df_gene <- norm_expr_df_gene[, CT.keep]
## AUC_REGION
avg_expr <- AUC_REGION_sub.m %>%
  t() %>%
  as.data.frame() %>%
  mutate(celltype = cell_metadata_sub$anno_wnn_v51) %>%
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
norm_expr_df_region <- norm_expr_df_region[, CT.keep]

# Get Gene Order
norm_expr_df_gene2 <- as.matrix(norm_expr_df_TF)
norm_expr_df_gene2 <- t(apply(norm_expr_df_gene2, 1, function(x) {
  (x - min(x)) / (max(x) - min(x))
}))
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
# df_color <- melt(tmp, varnames = c("Gene", "CellType"), value.name = "TF")
df_color <- as.data.frame(tmp, check.names = FALSE) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "CellType", values_to = "TF")
# df_size  <- melt(as.matrix(norm_expr_df_gene),  varnames = c("Gene", "CellType"), value.name = "Gene_base")
df_size <- as.data.frame(norm_expr_df_gene, check.names = FALSE) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "CellType", values_to = "Gene_base")
# df_alpha <- melt(as.matrix(norm_expr_df_region), varnames = c("Gene", "CellType"), value.name = "Region_base")
df_alpha <- as.data.frame(norm_expr_df_region, check.names = FALSE) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "CellType", values_to = "Region_base")
df_all <- purrr::reduce(list(df_color, df_size, df_alpha), dplyr::left_join, by = c("Gene", "CellType"))

# Generate Plots
## With TF Expression
Cairo::CairoPDF(paste0(FIG_DIR, "FL_selected_AUCell_heatmap.pdf"), 
                width = 13, height = 8, family = "Arial")
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
  scale_x_discrete(limits = rev(CT.keep)) +
  labs(y = "TF", x = "Cell Types", title = "TF Expr & AUCell Scores for HSPCs in PCW5") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = .5, colour = "black", size = 17)
  ) +
  coord_flip()
dev.off()


