###############################################################################
# Generate UMAP and Violin plots for customized region/gene set AUCell scores #
###############################################################################

##################
# Initialization #
##################

setwd("~/work/")
dir.create("./plots", showWarnings = F)
dir.create("./results", showWarnings = F)
source("./00.Initialization.R")

#############
# Load data #
#############

# Load qHSC AUCell region-based results
Region_based.h5 <- LoadH5Seurat("../Region_based_AUCell.h5seurat")

# Load Dev_M scRNA-seq object
FL.SeuratObj <- readRDS("~/local_data/FL_scrna_seurat.rds")

## Add AUCell data to Seurat Object metadata
cells_AUC.m <- GetAssayData(Region_based.h5)
cells_AUC.m <- cells_AUC.m[, colnames(FL.SeuratObj)]
cells_AUC.df <- as.data.frame(t(cells_AUC.m))
FL.SeuratObj <- AddMetaData(
  object = FL.SeuratObj,
  metadata = cells_AUC.df
)

#################################
# Generate Violinplots/boxplots #
#################################

AUC.idx <- colnames(cells_AUC.df)
gene.idx <- gsub(AUC.idx, pattern = "-region", replacement = "", fixed = TRUE) |>
  intersect(rownames(FL.SeuratObj))
ct_col <- "anno_wnn_v51"
interested_TF.idx <- c("GATA1", "GATA2", "GATA4", "GATA5")
interested_TF_AUC.idx <- paste0(interested_TF.idx, "-region")

# AUCell
for (i in interested_TF_AUC.idx) {
  idx_name <- i
  
  df <- FetchData(FL.SeuratObj, vars = c(ct_col, idx_name))
  colnames(df) <- c("CellType", "AUC")
  df <- df |> filter(is.finite(AUC))
  
  lev <- df |>
    group_by(CellType) |>
    summarize(stat = mean(AUC, na.rm = TRUE), .groups = "drop") |>
    arrange(desc(stat)) |>
    pull(CellType)
  
  df$CellType <- factor(df$CellType, levels = lev)
  
  p <- ggplot(df, aes(x = CellType, y = AUC, fill = CellType)) +
    geom_violin(scale = "width", trim = TRUE, color = "grey25") +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7, color = "black") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 1.6,
                 fill = "white", color = "black") +
    scale_fill_manual(values = COL) +
    labs(
      title = paste0("AUCell Scores by Cell Types: ", idx_name),
      x = "Cell Types", y = "AUCell Scores"
    ) +
    theme_classic() +
    theme(axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1), 
          axis.text = element_text(face = "bold"),
          panel.grid = element_blank(),
          plot.title = element_text(face = "bold", hjust = .5, colour = "black", size = 17)) +
    guides(fill = "none")
  
  save_plot(
    filename = paste0("./plots/",
                      idx_name, "_AUC_violin.png"),
    p, base_height = 6, base_width = 12, dpi = 300
  )
}

# Expression
for (i in interested_TF.idx) {
  idx_name <- i
  
  df <- FetchData(FL.SeuratObj, vars = c(ct_col, idx_name), layer = "data")
  colnames(df) <- c("CellType", "Expression")
  df <- df |> filter(is.finite(Expression))
  
  lev <- df |>
    group_by(CellType) |>
    summarize(stat = mean(Expression, na.rm = TRUE), .groups = "drop") |>
    arrange(desc(stat)) |>
    pull(CellType)
  
  df$CellType <- factor(df$CellType, levels = lev)
  
  p <- ggplot(df, aes(x = CellType, y = Expression, fill = CellType)) +
    geom_violin(scale = "width", trim = TRUE, color = "grey25") +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7, color = "black") +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 1.6,
                 fill = "white", color = "black") +
    scale_fill_manual(values = COL) +
    labs(
      title = paste0("Expression Level by Cell Types: ", idx_name),
      x = "Cell Types", y = "Expression"
    ) +
    theme_classic() +
    theme(axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1), 
          axis.text = element_text(face = "bold"),
          panel.grid = element_blank(),
          plot.title = element_text(face = "bold", hjust = .5, colour = "black", size = 17)) +
    guides(fill = "none")
  
  save_plot(
    filename = paste0("./plots/",
                      idx_name, "_Exp_violin.png"),
    p, base_height = 6, base_width = 12, dpi = 100
  )
}

###########################
# Generate WNN UMAP plots #
###########################

# For AUCell scores
for (i in interested_TF_AUC.idx) {
  idx_name <- i
  
  umap_df <- Embeddings(FL.SeuratObj, "wnn") %>% 
    as.data.frame() %>% 
    rownames_to_column("cell") %>%
    mutate(AUC = FL.SeuratObj@meta.data[cell, idx_name]) %>%
    arrange(AUC)
  
  p <- ggplot(umap_df, aes(x = wnnUMAP_1, y = wnnUMAP_2, color = AUC)) +
    geom_point(size = 0.5) +
    scale_color_viridis_c(option = "magma", name = idx_name) +
    theme_void() +
    theme(
      legend.position = "right",
      plot.title      = element_text(face = "bold", hjust = 0.5)
    ) +
    ggtitle(paste0(idx_name, "\nAUCell scores"))
  
  # p
  
  save_plot(filename = paste0("./plots/", idx_name, "_AUC_UMAP.png"),
            p, base_height = 10, base_width = 11, dpi =300)
}

# Expression
exp.m <- GetAssayData(FL.SeuratObj, assay = "RNA", layer = "data")
for (i in interested_TF.idx) {
  idx_name <- i
  
  umap_df2 <- Embeddings(FL.SeuratObj, "wnn") %>%
    as.data.frame() %>%
    rownames_to_column("cell") %>%
    mutate(Exp = as.vector(exp.m[idx_name ,])) %>%
    arrange(Exp)
  
  p <- ggplot(umap_df2, aes(x = wnnUMAP_1, y = wnnUMAP_2, color = Exp)) +
    geom_point(size = 0.5) +
    scale_color_viridis_c(option = "magma", name = idx_name) +
    theme_void() +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", hjust = 0.5)
    ) +
    ggtitle(paste0(idx_name, "\nExpression"))
  
  # p
  
  save_plot(filename = paste0("./plots/", idx_name, "_Exp_UMAP.png"),
            p, base_height = 10, base_width = 11, dpi=100)
}

