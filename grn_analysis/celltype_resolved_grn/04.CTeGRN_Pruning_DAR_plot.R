#############################
### Pruning CT-eGRN Plots ###
#############################

###-------------------Note-------------------###
###-------      Run with R-4.4.2      -------###
###------------------------------------------###

##################
# Initialization #
##################

setwd("~/local_data/proj/Dev_Multiome/04.regulome_R/01.SCENICplus/04.SCENICplus_CTeGRN/")
source("./00.Initialization.R")

library(AUCell)
library(ArchR)
library(caret)
library(pROC)
library(data.table)
library(future)
library(furrr)

# Setup Jaccard Index Similarity Function
jaccard_index <- function(A, B) {
  length(intersect(A, B)) / length(union(A, B))
}

prune_meta.l <- readRDS("./results/eRegulons_CT_Prune/CT_eGRN_DAR_Prune_Metadata.rds")
prune_meta.df <- lapply(prune_meta.l, function(x) x[[1]])
prune_meta.df <- bind_rows(prune_meta.df)
prune_meta_TF.all <- prune_meta.df %>% pull(TF) %>% unique

##################
# Generate Plots #
##################

##############################################
# Check CT-eGRN similarity across cell types #
##############################################

# Get pruned eRegulon
celltype_TF_prune <- list()
celltype_TG_prune <- list()
celltype_TR_prune <- list()
for (CT in CT_ORDER) {
  celltype_TF_prune[[CT]] <- prune_meta.l[[CT]]$CT_eRegulon_prune_meta %>%
    pull(TF) %>%
    unique()
  celltype_TG_prune[[CT]] <- prune_meta.l[[CT]]$CT_eRegulon_prune_meta %>%
    pull(Gene) %>%
    unique()
  celltype_TR_prune[[CT]] <- prune_meta.l[[CT]]$CT_eRegulon_prune_meta %>%
    pull(Region) %>%
    unique()
}
jaccard_matrix_TF <- outer(names(celltype_TF_prune), names(celltype_TF_prune),
                           Vectorize(function(x, y) jaccard_index(celltype_TF_prune[[x]], celltype_TF_prune[[y]])))
rownames(jaccard_matrix_TF) <- colnames(jaccard_matrix_TF) <- names(celltype_TF_prune)
jaccard_matrix_TG <- outer(names(celltype_TG_prune), names(celltype_TG_prune),
                           Vectorize(function(x, y) jaccard_index(celltype_TG_prune[[x]], celltype_TG_prune[[y]])))
rownames(jaccard_matrix_TG) <- colnames(jaccard_matrix_TG) <- names(celltype_TG_prune)
jaccard_matrix_TR <- outer(names(celltype_TR_prune), names(celltype_TR_prune),
                           Vectorize(function(x, y) jaccard_index(celltype_TR_prune[[x]], celltype_TR_prune[[y]])))
rownames(jaccard_matrix_TR) <- colnames(jaccard_matrix_TR) <- names(celltype_TR_prune)
## Plot
### Hierachical clustering
dist_mat <- as.dist(1 - jaccard_matrix_TR)
hc <- hclust(dist_mat, method = "complete")
celltype_order <- rownames(jaccard_matrix_TF)[hc$order]
### Construct long data.frame
jaccard_to_long <- function(jaccard_matrix, condition_name) {
  as.data.frame(as.table(jaccard_matrix)) %>%
    dplyr::rename(., CellType1 = Var1, CellType2 = Var2, Jaccard = Freq) %>%
    mutate(Condition = condition_name)
}
df_TF <- jaccard_to_long(jaccard_matrix_TF, "TF")
df_TG <- jaccard_to_long(jaccard_matrix_TG, "TG")
df_TR <- jaccard_to_long(jaccard_matrix_TR, "TR")
jaccard_df <- df_TF %>%
  dplyr::select(CellType1, CellType2, Jaccard_TF = Jaccard) %>%
  left_join(df_TG %>% dplyr::select(CellType1, CellType2, Jaccard_TG = Jaccard), by = c("CellType1", "CellType2")) %>%
  left_join(df_TR %>% dplyr::select(CellType1, CellType2, Jaccard_TR = Jaccard), by = c("CellType1", "CellType2"))
### Re-order cell type
jaccard_df$CellType1 <- factor(jaccard_df$CellType1, levels = celltype_order)
jaccard_df$CellType2 <- factor(jaccard_df$CellType2, levels = celltype_order)
# Only keep upper triangle
jaccard_df <- jaccard_df %>%
  filter(as.numeric(CellType1) >= as.numeric(CellType2))
### Plot
Cairo::CairoPDF(paste0(FIG_DIR, "04.CTeGRN_Pruning/CT_eGRN_DAR_Prune/CT_eGRN_DAR_Prune_JI_TR_heatmap.pdf"), 
                width = 15, height = 12)
ggplot(jaccard_df, aes(CellType1, CellType2)) +
  geom_tile(aes(fill = Jaccard_TF), na.rm = TRUE) +
  geom_point(aes(size = Jaccard_TG, color = Jaccard_TR), stroke = .5) +
  scale_fill_gradient(low = "white", high = "#E15759", na.value = NA) +
  scale_color_gradient(low = "white", high = "black") +
  labs(title = "Jaccard Similarity by Cell Type Specifc eGRNs\n(Based on Target Region Similarity)",
       fill = "Jaccard by TF (Heatmap Color)",
       size = "Jaccard by TG (Dot Size)",
       color = "Jaccard by TR (Dot Transparency)") +
  scale_y_discrete(position = "right") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.ticks.x = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(face = "bold"),
        panel.grid = element_blank(),
        legend.position = "left",
        plot.title = element_text(face = "bold", hjust = .5, colour = "black", size = 24))
dev.off()

#############################################
# Differential Target JI plots for every TF #
#############################################

# For all eRegulon, check JI across cell type
for (TF.interest in prune_meta_TF.all) {
  eRegulon_prune_TG <- list()
  eRegulon_prune_TR <- list()
  for (CT in CT_ORDER) {
    eRegulon_prune_TG[[CT]] <- prune_meta.l[[CT]]$CT_eRegulon_prune_meta %>%
      filter(TF == TF.interest) %>%
      pull(Gene) %>%
      unique()
    eRegulon_prune_TR[[CT]] <- prune_meta.l[[CT]]$CT_eRegulon_prune_meta %>%
      filter(TF == TF.interest) %>%
      pull(Region) %>%
      unique()
  }
  eRegulon_prune_TG <- Filter(length, eRegulon_prune_TG)
  eRegulon_prune_TR <- Filter(length, eRegulon_prune_TR)
  JI_eRegulon_TG.m <- outer(names(eRegulon_prune_TG), names(eRegulon_prune_TG),
                            Vectorize(function(x, y) jaccard_index(eRegulon_prune_TG[[x]], eRegulon_prune_TG[[y]])))
  rownames(JI_eRegulon_TG.m) <- colnames(JI_eRegulon_TG.m) <- names(eRegulon_prune_TG)
  JI_eRegulon_TR.m <- outer(names(eRegulon_prune_TR), names(eRegulon_prune_TR),
                            Vectorize(function(x, y) jaccard_index(eRegulon_prune_TR[[x]], eRegulon_prune_TR[[y]])))
  rownames(JI_eRegulon_TR.m) <- colnames(JI_eRegulon_TR.m) <- names(eRegulon_prune_TR)
  
  if (nrow(JI_eRegulon_TR.m) > 1) {
    ## Plot
    ### Hierachical clustering
    hc2 <- hclust(dist(JI_eRegulon_TR.m))
    # hc2 <- hclust(dist(JI_eRegulon_TG.m))
    CT_order <- rownames(JI_eRegulon_TR.m)[hc2$order]
    ### Construct long data.frame
    jaccard_to_long <- function(jaccard_matrix, condition_name) {
      as.data.frame(as.table(jaccard_matrix)) %>%
        dplyr::rename(., CellType1 = Var1, CellType2 = Var2, Jaccard = Freq) %>%
        mutate(Condition = condition_name)
    }
    df_eRegulonTG <- jaccard_to_long(JI_eRegulon_TG.m, "TG")
    df_eRegulonTR <- jaccard_to_long(JI_eRegulon_TR.m, "TR")
    jaccard_eRegulon_df <- df_eRegulonTG %>%
      dplyr::select(CellType1, CellType2, Jaccard_TG = Jaccard) %>%
      left_join(df_eRegulonTR %>% dplyr::select(CellType1, CellType2, Jaccard_TR = Jaccard), by = c("CellType1", "CellType2"))
    ### Re-order cell type
    jaccard_eRegulon_df$CellType1 <- factor(jaccard_eRegulon_df$CellType1, levels = CT_order)
    jaccard_eRegulon_df$CellType2 <- factor(jaccard_eRegulon_df$CellType2, levels = CT_order)
    # Only keep upper triangle
    jaccard_eRegulon_df <- jaccard_eRegulon_df %>%
      filter(as.numeric(CellType1) >= as.numeric(CellType2))
    ### Plot
    p <- ggplot(jaccard_eRegulon_df, aes(CellType1, CellType2)) +
      # geom_tile(aes(fill = Jaccard_TF), na.rm = TRUE) +
      geom_point(aes(size = Jaccard_TG, color = Jaccard_TR), stroke = .5) +
      # scale_fill_gradient(low = "white", high = "#E15759", na.value = NA) +
      scale_color_gradient(low = "#E0FFFF", high = "red", limits = c(0,1)) +
      scale_size_continuous(limits = c(0,1), range = c(1, 6)) +
      theme_classic() +
      labs(title = paste0("Jaccard Similarity on ", TF.interest, " eRegulon"),
           fill = "Jaccard by TF (Heatmap)",
           size = "Jaccard by TG (Dot Size)",
           color = "Jaccard by TR (Dot Color)") +
      scale_y_discrete(position = "right") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(face = "bold", hjust = .5, colour = "black", size = 17),
            legend.position = "left",
            axis.title = element_blank())
    save_plot(filename = paste0("./plots/04.CTeGRN_Pruning/CT_eGRN_DAR_Prune_new/TF_JI_Heatmaps/", TF.interest, "_JI.pdf"),
              p, base_height = 10, base_width = 10)
  }
}

##############################################
# Compute and plot differential target ratio #
##############################################

# Compute differential target ratio
DT.l <- list()
for (TF.interest in prune_meta_TF.all) {
  ## Compute JI
  eRegulon_prune_TG <- list()
  for (CT in CT_ORDER) {
    eRegulon_prune_TG[[CT]] <- prune_meta.l[[CT]]$CT_eRegulon_prune_meta %>%
      filter(TF == TF.interest) %>%
      pull(Gene) %>%
      unique()
  }
  eRegulon_prune_TG <- Filter(length, eRegulon_prune_TG)
  JI_eRegulon_TG.m <- outer(names(eRegulon_prune_TG), names(eRegulon_prune_TG),
                            Vectorize(function(x, y) jaccard_index(eRegulon_prune_TG[[x]], eRegulon_prune_TG[[y]])))
  rownames(JI_eRegulon_TG.m) <- colnames(JI_eRegulon_TG.m) <- names(eRegulon_prune_TG)
  if (nrow(JI_eRegulon_TG.m) > 1) {
    ## Separate differential targets by walktrap algorithm
    library(igraph)
    g <- graph_from_adjacency_matrix(
      adjmatrix = JI_eRegulon_TG.m,
      mode      = "undirected",
      weighted  = TRUE,
      diag      = FALSE
    )
    ### Run walktap
    wc <- cluster_walktrap(g, weights = E(g)$weight, steps = 4)
    n_comm <- length(wc)
    membership <- membership(wc)
    DT.l[[TF.interest]] <- list(TF.interest = TF.interest,
                                n_CT = nrow(JI_eRegulon_TG.m), 
                                n_comm = n_comm, 
                                membership = membership)
  }else {
    membership <- 1
    names(membership) <- rownames(JI_eRegulon_TG.m)
    DT.l[[TF.interest]] <- list(TF.interest = TF.interest,
                                n_CT = 1, 
                                n_comm = 1, 
                                membership = membership)
  }
}

df_points <- map_df(
  DT.l,
  ~ tibble(
    TF = .x[[1]],
    x = .x[[2]],
    y = .x[[3]]
  )
)

df_plot <- df_points %>%
  mutate(y = as.factor(as.character(y)))

df_jit <- df_plot %>%
  mutate(
    # map y‐factor → 1,2,3,… 
    grp_num = as.numeric(y),
    # add a single random jitter to each row
    grp_jit = grp_num + runif(n(), -0.2, 0.2)
  )

p <- ggplot(df_jit, aes(x = grp_jit, y = x)) +
  # points using that jit‐x
  geom_point(aes(color = y), size = 1.5, alpha = 0.7, position = position_identity()) +
  # text labels using the same jit‐x
  geom_text_repel(
    aes(label = TF, color = y),
    position      = position_identity(),
    size          = 3,
    box.padding   = 0.2, 
    force         = 2,
    point.padding = 0.2,
    segment.size  = 0.3,
    max.overlaps  = Inf,
    seed          = 42
  ) +
  # now draw a continuous x‐axis but label positions 1,2,3 with your factor levels
  scale_x_continuous(
    breaks = seq_along(levels(df_jit$y)),
    labels = levels(df_jit$y),
    expand = c(0.2, 0.2)
  ) +
  scale_color_manual(values = COL4) +
  guides(color = "none") +
  labs(
    title = "Differential Target Analysis of CT-eRegulons",
    x     = "Number of Modules Deteected with Jaccard Similarity Network",
    y     = "Number of cell types containing this TF in their CT-eRegulons"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title  = element_text(face = "bold", hjust = 0.5),
    axis.title  = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold")
  )
p
save_plot(filename = paste0("./plots/04.CTeGRN_Pruning/CT_eGRN_DAR_Prune_new/Differential_targets_Walktrap_Summary.pdf"),
          p, base_height = 10, base_width = 12)

# New version of diff target sum plot
df_plot <- df_points %>%
  transmute(
    TF,
    n_CT   = x,  # Number of cell types
    n_comm = y   # Number of modules (communities)
  )

x_cut <- median(df_plot$n_comm, na.rm = TRUE)
y_cut <- median(df_plot$n_CT,   na.rm = TRUE)

TF_highlight <- c(
  "RUNX1", "MYB", "HMGA2", "GATA2", "MECOM",
  "LIN28B", "ETV6", "SMARCC1", "TBL1XR1"
)

df_plot <- df_plot %>%
  mutate(
    n_comm_f = factor(n_comm), 
    TF_lab   = if_else(TF %in% TF_highlight, TF, NA_character_)
  )

p <- ggplot(df_plot, aes(x = n_comm, y = n_CT)) +
  stat_density_2d_filled(
    geom        = "polygon",
    contour_var = "ndensity",
    bins        = 30,
    alpha       = 0.4,
    show.legend = FALSE
  ) +
  geom_vline(xintercept = x_cut, linetype = "dashed", linewidth = 0.4, colour = "grey30") +
  geom_hline(yintercept = y_cut, linetype = "dashed", linewidth = 0.4, colour = "grey30") +
  geom_point(
    aes(colour = n_comm_f, size = n_CT),
    alpha = 0.8
  ) +
  geom_text_repel(
    data         = subset(df_plot, !is.na(TF_lab)),
    aes(label    = TF_lab, colour = n_comm_f),
    size         = 3,
    box.padding  = 0.25,
    point.padding= 0.2,
    segment.size = 0.3,
    max.overlaps = Inf,
    seed         = 42
  ) +
  scale_x_continuous(
    breaks = pretty(range(df_plot$n_comm)),
    minor_breaks = NULL
  ) +
  scale_y_continuous(
    breaks = pretty(range(df_plot$n_CT)),
    minor_breaks = NULL
  ) +
  scale_colour_manual(
    values = COL4,
    name   = "Number of modules"
  ) +
  scale_size_continuous(
    range = c(1, 10),
    guide = "none"
  ) +
  labs(
    title = "Differential Target Analysis of CT-eRegulons",
    x     = "Number of modules detected with Jaccard similarity network",
    y     = "Number of cell types containing this TF in their CT-eRegulons"
  ) +
  theme_classic() +
  theme(
    plot.title  = element_text(face = "bold", hjust = 0.5),
    axis.title  = element_text(face = "bold"),
    axis.text   = element_text(face = "plain")
  )
p

##########################################################
# Compute AUCell scores from interested pruned eRegulons #
##########################################################

# Load packages
library(AUCell)
library(GSEABase)

# Load snRNA-seq data with WNN UMAP
FL.SeuratObj <- readRDS("~/local_data/proj/Dev_Multiome/data/FL_scrna_seurat_20250401.rds")
FL_wnn.df <- read.csv("/work/DevM_analysis/01.annotation/10.integration_joint_clean/data/FL_wnn_umap.csv", 
                      row.names = 1)
rownames(FL_wnn.df) <- gsub(pattern = "_", replacement = "#", 
                            x = rownames(FL_wnn.df), fixed = TRUE)
FL_wnn.df <- FL_wnn.df[colnames(FL.SeuratObj) ,]
wnn_reduction <- CreateDimReducObject(embeddings = as.matrix(FL_wnn.df),
                                      key = "wnn_",
                                      assay = DefaultAssay(FL.SeuratObj))
FL.SeuratObj[["wnn"]] <- wnn_reduction

# ## Set EXP data & prepare ranking
# cell_meata <- FL.SeuratObj@meta.data
# exp.m <- GetAssayData(FL.SeuratObj, assay = "RNA")
# cells_rankings <- AUCell_buildRankings(exp.m, plotStats=FALSE)
# saveRDS(cells_rankings, file = "./AUCell_ranking.rds")
cells_rankings <- readRDS("./results/eRegulons_AUCell/GeneRank_for_AUCell.rds")
cells_rankings2 <- cells_rankings[, colnames(FL.SeuratObj)]

# Set up eRegulon TG sets
TF.interest <- "RUNX1"
## Compute JI
eRegulon_prune_TG <- list()
for (CT in CT_ORDER) {
  eRegulon_prune_TG[[CT]] <- prune_meta.l[[CT]]$CT_eRegulon_prune_meta %>%
    filter(TF == TF.interest) %>%
    pull(Gene) %>%
    unique()
}
eRegulon_prune_TG <- Filter(length, eRegulon_prune_TG)
JI_eRegulon_TG.m <- outer(names(eRegulon_prune_TG), names(eRegulon_prune_TG),
                          Vectorize(function(x, y) jaccard_index(eRegulon_prune_TG[[x]], eRegulon_prune_TG[[y]])))
rownames(JI_eRegulon_TG.m) <- colnames(JI_eRegulon_TG.m) <- names(eRegulon_prune_TG)
## Separate differential targets by walktrap algorithm
library(igraph)
g <- graph_from_adjacency_matrix(
  adjmatrix = JI_eRegulon_TG.m,
  mode      = "undirected",
  weighted  = TRUE,
  diag      = FALSE
)
### Run walktap
wc <- cluster_walktrap(g, weights = E(g)$weight, steps = 4)
n_comm <- length(wc)
membership <- membership(wc)
## Define gene set by memberships
geneSets.l <- list()
for (i in 1:n_comm) {
  CT.idx <- names(membership)[which(membership == i)]
  geneSets.tmp <- eRegulon_prune_TG[[CT.idx[1]]]
  for (j in CT.idx) {
    geneSets.tmp <- intersect(geneSets.tmp, eRegulon_prune_TG[[j]])
  }
  geneSets.l[[i]] <- geneSets.tmp
}
geneSets.l2 <- list()
for (i in 1:n_comm) {
  geneSets.tmp <- geneSets.l[[i]]
  for (j in setdiff(1:n_comm, i)) {
    geneSets.tmp <- setdiff(geneSets.tmp, geneSets.l[[j]])
  }
  geneSets.l2[[i]] <- GeneSet(geneSets.tmp, setName=as.character(paste0(TF.interest, "_", i)))
}
## Run AUCell
geneSets <- GeneSetCollection(geneSets.l2)
names(geneSets)
### Compute gene based AUCell score
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
saveRDS(cells_AUC, paste0("./results/eRegulons_CT_Prune/CT_eGRN_DAR_Prune_", TF.interest, "_AUCell.rds"))

# Add new AUCell data to Seurat Object
cells_AUC.m <- assay(cells_AUC)
rownames(cells_AUC.m) <- paste0(TF.interest, "_", 1:n_comm)
cells_AUC.m <- cells_AUC.m[, colnames(FL.SeuratObj)]
cells_AUC.df <- as.data.frame(t(cells_AUC.m))
FL.SeuratObj <- AddMetaData(
  object = FL.SeuratObj,
  metadata = cells_AUC.df
)
## Featureplot
idx <- paste0(TF.interest, "_", 2)
umap_df <- Embeddings(FL.SeuratObj, "wnn") %>% 
  as.data.frame() %>% 
  rownames_to_column("cell") %>%
  mutate(AUC = FL.SeuratObj@meta.data[cell, idx]) %>%
  arrange(AUC)

p <- ggplot(umap_df, aes(x = wnn_1, y = wnn_2, color = AUC)) +
  geom_point(size = 0.5) +
  scale_color_viridis_c(option = "magma", name = idx) +
  theme_void() +
  theme(
    legend.position = "right",
    plot.title      = element_text(face = "bold", hjust = 0.5)
  ) +
  ggtitle(paste0("UMAP colored by ", idx))

save_plot(filename = paste0("./plots/04.CTeGRN_Pruning/CT_eGRN_DAR_Prune/Differential_targets/", idx, "_UMAP.pdf"),
          p, base_height = 10, base_width = 10)

