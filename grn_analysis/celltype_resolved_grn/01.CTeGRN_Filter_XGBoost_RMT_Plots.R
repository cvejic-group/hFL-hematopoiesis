#############################################
### Check eGRN Filtering Results and Plot ###
#############################################

###-------------------Note-------------------###
###-------      Run with R-4.4.2      -------###
###------------------------------------------###

##################
# Initialization #
##################

setwd("~/work/")
source("./00.Initialization.R")

# Load additional packages
library(AUCell)
library(ArchR)

###################################################
### Load Cell Type specific eGRN by XGBoost RMT ###
###################################################

# Load CT-eRegulons info
CT_eReg_XGBoost_RMT.df <- readRDS(paste0(RES_DIR, "eRegulons_CT_Filter/XGBoost_RMT_TFRes.rds"))

# Form CT-eRegulons list
CT_eRegulon.l <- list()
CT_eReg.df <- CT_eReg_XGBoost_RMT.df
for (CT in CT_ORDER) {
  ### Add MAZ for Transient population
  if (CT %in% c("MEMP-t", "MastP-t", "MDP", "LMPP", "LP", "Cycling-LP")) {
    CT_eRegulon.l[[CT]] <- union(CT_eReg.df[[CT]]$confident_features, c("MAZ", "CHURC1", "CEBPZ"))
  } else {
    CT_eRegulon.l[[CT]] <- CT_eReg.df[[CT]]$confident_features 
  }
}

###########################
### Generate SHAP plots ###
###########################

# Setup function
plot_SHAP_vs_null <- function(shap_df, top_n = 50, title = "SHAP Rank Plot", label_all = FALSE) {
  shap_df <- shap_df[order(-shap_df$SHAP_ratio), ]
  shap_df$Rank <- seq_len(nrow(shap_df))
  
  top_df <- shap_df[1:top_n, ]
  top_df <- top_df[is.finite(top_df$SHAP) & is.finite(top_df$SHAP_random), ]
  
  # Define color group
  top_df$group <- "Excluded"
  top_df$group[top_df$CT_eReg_possible] <- "Possible"
  top_df$group[top_df$CT_eReg_confident] <- "Confident"
  top_df$group <- factor(top_df$group, levels = c("Confident", "Possible", "Excluded"))
  
  # Label
  if (label_all) {
    idx <- top_df$CT_eReg_confident + top_df$CT_eReg_possible > 0
    label_df <- top_df[idx ,]
  } else {
    idx <- top_df$CT_eReg_confident
    label_df <- top_df[idx ,]
  }
  
  top_df$SHAP_ratio_plot <- log1p(top_df$SHAP_ratio)
  top_df$SHAP_ratio_plot <- (top_df$SHAP_ratio_plot - min(top_df$SHAP_ratio_plot)) / (max(top_df$SHAP_ratio_plot) - min(top_df$SHAP_ratio_plot))
  
  p <- ggplot(top_df, aes(x = Rank)) +
    geom_point(aes(y = SHAP, color = group, size = SHAP_ratio_plot)) +
    geom_point(aes(y = SHAP_random), color = "darkgrey", shape = 1, size = 4) +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(y = SHAP, label = Feature),
      size = 5, max.overlaps = Inf, box.padding = 0.5,
      segment.size = 0.2, min.segment.length = 0.8
    ) +
    scale_color_manual(values = c("Confident" = "red", "Possible" = "orange", "Excluded" = "gray")) +
    labs(x = "TF (ranked by SHAP_ratio)", y = "SHAP score", title = title, color = "Category") +
    guides(size = "none") + 
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.text = element_text(face = "bold"),
          panel.grid = element_blank(),
          plot.title = element_text(face = "bold", hjust = .5, colour = "black", size = 17))
  return(p)
}

# Plot out
for (CT in CT_ORDER) {
  shap_df <- CT_eReg_XGBoost_RMT.df[[CT]]$shap
  p <- plot_SHAP_vs_null(shap_df, top_n = 148, 
                         title = paste0("SHAP Rank of ", CT), 
                         label_all = TRUE)
  
  ggsave(
    filename = paste0(FIG_DIR, "02.CTeGRN_Filter_XGBoost_RMT/SHAP_rank_plots/", 
                      CT, "_SHAP_Rank_plot.pdf"),
    plot = p,
    width = 16,
    height = 12,
    device = cairo_pdf,
    family = "Arial"
  )
}

##################################################
### Check the selected CT-eRegulons on heatmap ###
##################################################

# Load data
TOPeReg_TF_AUCell_mat <- readRDS(paste0(RES_DIR, "eRegulons_AUCell/TOPeReg_TF_AUCell_mat.rds"))
TF_EXPR.m <- TOPeReg_TF_AUCell_mat$TF_EXPR.m
AUC_GENE.m <- TOPeReg_TF_AUCell_mat$AUC_GENE.m
AUC_REGION.m <- TOPeReg_TF_AUCell_mat$AUC_REGION.m

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

# Prepare plot data
df_color <- melt(tmp, varnames = c("Gene", "CellType"), value.name = "TF")
df_size  <- melt(as.matrix(norm_expr_df_gene),  varnames = c("Gene", "CellType"), value.name = "Gene_base")
df_size$Gene <- df_color$Gene
df_alpha <- melt(as.matrix(norm_expr_df_region), varnames = c("Gene", "CellType"), value.name = "Region_base")
df_alpha$Gene <- df_color$Gene
df_all <- purrr::reduce(list(df_color, df_size, df_alpha), dplyr::left_join, by = c("Gene", "CellType"))
## Fetch CT-eRegulon
marker_df <- do.call(rbind, lapply(CT_ORDER, function(ct){
  data.frame(
    CellType = ct,
    Gene = CT_eRegulon.l[[ct]],
    stringsAsFactors = FALSE
  )
}))
marker_df$CellType <- factor(marker_df$CellType,
                             levels = CT_ORDER)
marker_df$Gene <- factor(marker_df$Gene,
                         levels = gene_order)
Cairo::CairoPDF(paste0(FIG_DIR, "02.CTeGRN_Filter_XGBoost_RMT/CT_eREG_Mark_heatmap.pdf"), 
                width = 48, height = 12, family = "Arial")
ggplot(df_all, aes(x = CellType, y = Gene)) +
  geom_tile(aes(fill = TF), color = "white") +
  geom_point(aes(size = Gene_base, alpha = Region_base), color = "black") +
  geom_tile(
    data  = marker_df,
    aes(x = CellType, y = Gene),
    fill  = NA,
    color = "black",
    size  = 0.8
  ) +
  scale_fill_gradientn(
    colours = c("#3a86ff", "#48cae4", "white", "#e9c46a", "#d62828"),
    values  = scales::rescale(c(-2, -1, 0, 1, 2)),
    name    = "TF expression"
  ) +
  scale_size(range = c(1, 5), name = "Gene-based AUCell scores") +
  scale_alpha(range = c(0.1, 1), name = "Region-based AUCell scores") +
  scale_y_discrete(limits = gene_order) +
  scale_x_discrete(limits = rev(CT_ORDER)) +
  labs(y = "TF", x = "Cell Types", title = "Mark CT-eRegulons for All Cell Types") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = .5, colour = "black", size = 17)
  ) +
  coord_flip()
dev.off()

######################################################
### Check Cell Type Similarity using Jaccard Index ###
######################################################

library(ggdendro)

# Setup Jaccard Index Similarity Function
jaccard_index <- function(A, B) {
  length(intersect(A, B)) / length(union(A, B))
}

# Compute JI
jaccard_matrix <- outer(names(CT_eRegulon.l), names(CT_eRegulon.l),
                        Vectorize(function(x, y) jaccard_index(CT_eRegulon.l[[x]], CT_eRegulon.l[[y]])))
rownames(jaccard_matrix) <- colnames(jaccard_matrix) <- names(CT_eRegulon.l)
## Clip values
jaccard_matrix.m <- jaccard_matrix
# jaccard_matrix.m[which(jaccard_matrix.m <= 0.15)] <- 0

# Hierarchical clustering
dist_mat <- as.dist(1 - jaccard_matrix.m)
hc <- hclust(dist_mat, method = "complete")
ordered_ct <- hc$labels[hc$order]
ordered_mat <- jaccard_matrix.m[ordered_ct, ordered_ct]
diag(ordered_mat) <- 0.7
## Cut tree
k <- 7
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
### Add rectangle for transient populations
# cluster_rects <- rbind(cluster_rects, c(8, 5.5, 9.5, 20.5, 22.5))
# cluster_rects <- rbind(cluster_rects, c(9, 11.5, 13.5, 20.5, 22.5))
# cluster_rects <- rbind(cluster_rects, c(10, 14.5, 16.5, 23.5, 25.5))
# cluster_rects <- rbind(cluster_rects, c(11, 14.5, 16.5, 27.5, 31.5))
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
  geom_rect(
    data = cluster_rects,
    inherit.aes = FALSE,
    aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax,
        colour = rect_color),
    fill = NA, 
    colour = cluster_rects$rect_color,
    linewidth = 1,
    show.legend = TRUE
  ) +
  scale_x_discrete(expand = c(0, 0), limits = ordered_ct) +
  scale_y_discrete(expand = c(0, 0), limits = rev(ordered_ct)) +
  # SCAV_COLS=colorRampPalette(c("#F9F6F4", "#F1DEDA", "#C684AA", "#6B4E80", "#4D4665"))(100)
  scale_fill_gradientn(
    # colours = c("#253494", "#1F78B4", "white", "#FF7F00", "#8B0000"),
    colours = c("#F9F6F4", "#F1DEDA", "#C684AA", "#6B4E80", "#4D4665"),
    # na.value = "lightgray",
    name = "Jaccard Index"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(face = "bold"),
        axis.title   = element_blank(),
        panel.grid = element_blank(),
        plot.margin  = margin(t = -2, b = 5, l = 5, r = 5),
        plot.title = element_text(face = "bold", hjust = .5, colour = "black", size = 17))

# Combine dendrogram and heatmap
Cairo::CairoPDF(paste0(FIG_DIR, "02.CTeGRN_Filter_XGBoost_RMT/CT_eGRN_Filter_JI_heatmap_woRED2.pdf"), 
                width = 12, height = 11, family = "Arial")
plot_grid(
  p_dendro, p_heat,
  ncol = 1,
  rel_heights = c(0.07, 1),
  align = "v"
)
dev.off()

###################################
### Save Filtered eGRN to table ###
###################################

# Get CT-eGRN Filtered
eGRN_filter_meta_df.l <- list()
for(CT in CT_ORDER){
  eGRN_filter_meta_df.l[[CT]] <- BA_eRegulon_meta.df %>%
    dplyr::filter(TF %in% CT_eRegulon.l[[CT]])
  eGRN_filter_meta_df.l[[CT]]$CT_eGRN_name <- paste0(CT, 
                                                     "_", 
                                                     eGRN_filter_meta_df.l[[CT]]$TF)
}
saveRDS(eGRN_filter_meta_df.l, 
        paste0(RES_DIR, "eRegulons_CT_Filter/CT_eGRN_Filter_Metadata.rds"))
### Save to xlsx
library(openxlsx)
wb <- createWorkbook()
for(CT in CT_ORDER){
  addWorksheet(wb, CT)
  writeData(wb, CT, eGRN_filter_meta_df.l[[CT]])
}
saveWorkbook(wb, paste0(RES_DIR, "eRegulons_CT_Filter/CT_eGRN_Filter_Metadata.xlsx"), 
             overwrite = TRUE)



