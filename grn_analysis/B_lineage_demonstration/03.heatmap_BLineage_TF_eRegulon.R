####################################################
### Filtered CT-eRegulon B Lineage TF similarity ###
####################################################

###-------------------Note-------------------###
###-------      Run with R-4.4.2      -------###
###------------------------------------------###

##################
# Initialization #
##################

setwd("~/")
source("./00.Initialization.R")
dir.create(paste0(FIG_DIR, "03.heatmap_BLineage_TF_eRegulon"))

# Load additional packages
library(AUCell)
library(ArchR)

###################################################
### Load Cell Type specific eGRN by XGBoost RMT ###
###################################################

XGBoost_RMT.o <- read_rds("./data/XGBoost_RMT_TFRes_wRankedFeature.rds")

CT_eRegulon.l <- list()
for (CT in CT_ORDER) {
  CT_eRegulon.l[[CT]] <- XGBoost_RMT.o[[CT]]$confident_features
}

CT.v <- c("HSC", "PreProB", 
          "ProB-1",
          "ProB-2",
          "Large-PreB",
          "Small-PreB", 
          "IM-B")

BLineage_eRegulon.l <- CT_eRegulon.l[CT.v]

######################################################
### Check Cell Type Similarity using Jaccard Index ###
######################################################

library(ggdendro)

# Setup Jaccard Index Similarity Function
jaccard_index <- function(A, B) {
  length(intersect(A, B)) / length(union(A, B))
}

# Compute JI
jaccard_matrix <- outer(names(BLineage_eRegulon.l), names(BLineage_eRegulon.l),
                        Vectorize(function(x, y) jaccard_index(BLineage_eRegulon.l[[x]], BLineage_eRegulon.l[[y]])))
rownames(jaccard_matrix) <- colnames(jaccard_matrix) <- names(BLineage_eRegulon.l)
## Clip values
jaccard_matrix.m <- jaccard_matrix
# jaccard_matrix.m[which(jaccard_matrix.m <= 0.15)] <- 0

# Hierarchical clustering
dist_mat <- as.dist(1 - jaccard_matrix.m)
hc <- hclust(dist_mat, method = "complete")
## Override the cluster order
ordered_ct <- CT.v
# ordered_ct <- hc$labels[hc$order]
ordered_mat <- jaccard_matrix.m[ordered_ct, ordered_ct]
diag(ordered_mat) <- 0.7
## Cut tree
k <- 2
clusters <- cutree(hc, k = k)
cluster_rects <- tibble(
  CellType = names(clusters),
  cluster  = clusters
) %>%
  mutate(CellType = factor(CellType, levels = ordered_ct)) %>%
  group_by(cluster) %>%
  summarise(
    start = min(as.integer(CellType)) - 0.5,
    end   = max(as.integer(CellType)) + 0.5,
    .groups = "drop"
  ) %>%
  mutate(
    ## for the heatmap y-axis we reversed the order
    ymin = length(ordered_ct) - end + 1,
    ymax = length(ordered_ct) - start + 1
  )
### Set color
cluster_rects <- cluster_rects %>%
  mutate(
    rect_color = ifelse(row_number() <= 7, "black", "red")
  )

# Prepare plot data
ddata <- dendro_data(as.dendrogram(hc), type = "rectangle")
df <- as.data.frame(ordered_mat) %>%
  rownames_to_column("CellType1") %>%
  pivot_longer(-CellType1, names_to = "CellType2", values_to = "Jaccard") %>%
  mutate(
    CellType1 = factor(CellType1, levels = rev(ordered_ct)),
    CellType2 = factor(CellType2, levels = ordered_ct)
  )

# Plot dendrogram
p_dendro <- ggplot(segment(ddata), 
                   aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_segment() +
  scale_x_continuous(
    limits = c(0.5, length(ordered_ct) + 0.5),
    expand = c(0, 0)
  ) +
  theme_void() +
  labs(title = "Cell Type Similarity by TFs") +
  theme(plot.margin  = margin(t = -2, b = 5, l = 5, r = 5),
        plot.title = element_text(face = "bold", hjust = .5, colour = "black", size = 17))

# Plot heatmap with cluster rectangles
p_heat <- ggplot(df, aes(x = CellType2, y = CellType1, fill = Jaccard)) +
  geom_tile(color = "white") +
  
  scale_x_discrete(expand = c(0, 0), limits = ordered_ct) +
  scale_y_discrete(expand = c(0, 0), limits = rev(ordered_ct)) +
  scale_fill_gradientn(
    colours = c("#F9F6F4", "#F1DEDA", "#C684AA", "#6B4E80", "#4D4665"),
    name = "Jaccard Index"
  ) +
  labs(title = "Cell Type Similarity by TFs") + 
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(face = "bold"),
        axis.title   = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = .5, colour = "black", size = 17))

# plot
Cairo::CairoPDF(paste0(FIG_DIR, "03.heatmap_BLineage_TF_eRegulon/BLineage_TF_JI_heatmap2.pdf"), 
                width = 8, height = 7, family = "Arial")
p_heat
dev.off()

##################################################
### Check the selected CT-eRegulons on heatmap ###
##################################################

# Load data
cell_metadata <- read_rds("./data/FL_cell_metadata.rds")
TOPeReg_TF_AUCell_mat <- readRDS("./data/TOPeReg_TF_AUCell_mat.rds")
TF_EXPR.m <- TOPeReg_TF_AUCell_mat$TF_EXPR.m
AUC_GENE.m <- TOPeReg_TF_AUCell_mat$AUC_GENE.m
AUC_REGION.m <- TOPeReg_TF_AUCell_mat$AUC_REGION.m

# TF.v <- unlist(CT_eRegulon.l) |> unique()

TF.v <- c("HLF", "MYCN", "HMGA2", "MEIS1", "MECOM", "NFATC2",
          "FOXO1", "TCF3", "EBF1", "PAX5")

# Compute AVG
## EXPR
avg_expr <- TF_EXPR.m[TF.v, ] %>%
  t() %>%
  as.data.frame() %>%
  mutate(celltype = cell_metadata$anno_wnn_v51) %>%
  group_by(celltype) %>%
  filter(celltype %in% CT.v) %>%
  summarise(across(everything(), mean)) %>%
  column_to_rownames(var = "celltype") %>%
  t()
zscore_expr <- t(apply(avg_expr, 1, function(x) {
  (x - mean(x)) / sd(x)
}))
norm_expr_df_TF <- as.data.frame(zscore_expr)
norm_expr_df_TF <- norm_expr_df_TF[, CT.v]
norm_expr_df_TF2 <- t(apply(zscore_expr, 1, function(x) {
  (x - min(x)) / (max(x) - min(x))
}))
norm_expr_df_TF2 <- as.data.frame(norm_expr_df_TF2)
norm_expr_df_TF2 <- norm_expr_df_TF2[, CT.v]
## AUC_GENE
avg_expr <- AUC_GENE.m[TF.v, ] %>%
  t() %>%
  as.data.frame() %>%
  mutate(celltype = cell_metadata$anno_wnn_v51) %>%
  group_by(celltype) %>%
  filter(celltype %in% CT.v) %>%
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
norm_expr_df_gene <- norm_expr_df_gene[, CT.v]
## AUC_REGION
avg_expr <- AUC_REGION.m[TF.v ,] %>%
  t() %>%
  as.data.frame() %>%
  mutate(celltype = cell_metadata$anno_wnn_v51) %>%
  group_by(celltype) %>%
  filter(celltype %in% CT.v) %>%
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
norm_expr_df_region <- norm_expr_df_region[, CT.v]

# Get Gene Order
# norm_expr_df_gene2 <- as.matrix(norm_expr_df_gene)
norm_expr_df_gene2 <- as.matrix(norm_expr_df_TF2)
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

# Prepare plot data
df_color <- reshape2::melt(tmp, varnames = c("Gene", "CellType"), value.name = "TF")
df_size  <- reshape2::melt(as.matrix(norm_expr_df_gene),  varnames = c("Gene", "CellType"), value.name = "Gene_base")
df_size$Gene <- df_color$Gene
df_alpha <- reshape2::melt(as.matrix(norm_expr_df_region), varnames = c("Gene", "CellType"), value.name = "Region_base")
df_alpha$Gene <- df_color$Gene
df_all <- purrr::reduce(list(df_color, df_size, df_alpha), dplyr::left_join, by = c("Gene", "CellType"))
## Fetch CT-eRegulon
marker_df <- do.call(rbind, lapply(CT.v, function(ct){
  data.frame(
    CellType = ct,
    Gene = CT_eRegulon.l[[ct]],
    stringsAsFactors = FALSE
  )
}))
marker_df$CellType <- factor(marker_df$CellType,
                             levels = CT.v)
marker_df$Gene <- factor(marker_df$Gene,
                         levels = gene_order)
Cairo::CairoPDF("./plots/03.heatmap_BLineage_TF_eRegulon/BLineage_Markers_heatmap2.pdf", 
                width = 8, height = 6, family = "Arial")
ggplot(df_all, aes(x = CellType, y = Gene)) +
  geom_tile(aes(fill = TF), color = "white") +
  geom_point(aes(size = Gene_base, alpha = Region_base), color = "black") +
  geom_tile(
    data  = marker_df,
    aes(x = CellType, y = Gene),
    fill  = NA,
    color = "black",
    linewidth  = 0.8
  ) +
  scale_fill_gradientn(
    colours = c("#3a86ff", "#48cae4", "white", "#e9c46a", "#d62828"),
    values  = scales::rescale(c(-2, -1, 0, 1, 2)),
    name    = "TF expression"
  ) +
  scale_size(range = c(1, 5), name = "Gene-based AUCell scores") +
  scale_alpha(range = c(0.1, 1), name = "Region-based AUCell scores") +
  scale_y_discrete(limits = gene_order) +
  scale_x_discrete(limits = rev(CT.v)) +
  labs(y = "", x = "", title = "") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = .5, colour = "black", size = 17),
    legend.position = "bottom",
    legend.box = "vertical"
  ) +
  coord_flip()
dev.off()

