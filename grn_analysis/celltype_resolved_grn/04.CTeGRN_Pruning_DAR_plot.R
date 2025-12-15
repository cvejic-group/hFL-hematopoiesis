#############################
### Pruning CT-eGRN Plots ###
#############################

###-------------------Note-------------------###
###-------      Run with R-4.4.2      -------###
###------------------------------------------###

##################
# Initialization #
##################

setwd("~/work/")
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

prune_meta.l <- readRDS("CT_eGRN_DAR_Prune_Metadata.rds")
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
Cairo::CairoPDF(paste0(FIG_DIR, "CT_eGRN_DAR_Prune_JI_TR_heatmap.pdf"), 
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
    save_plot(filename = paste0("TF_JI_Heatmaps/", TF.interest, "_JI.pdf"),
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
    grp_num = as.numeric(y),
    grp_jit = grp_num + runif(n(), -0.2, 0.2)
  )

p <- ggplot(df_jit, aes(x = grp_jit, y = x)) +
  geom_point(aes(color = y), size = 1.5, alpha = 0.7, position = position_identity()) +
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
save_plot(filename = paste0("Differential_targets_Walktrap_Summary.pdf"),
          p, base_height = 10, base_width = 12)
