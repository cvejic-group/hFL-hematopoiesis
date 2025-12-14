#############################
### CT-eGRN Network Plots ###
#############################

###-------------------Note-------------------###
###-------      Run with R-4.4.2      -------###
###------------------------------------------###

##################
# Initialization #
##################

setwd("~/local_data/proj/Dev_Multiome/04.regulome_R/06.SCENICplus_CTeGRN/")
source("./00.Initialization.R")

library(purrr)
library(tibble)
library(igraph)
library(readr)

# Load CT-eGRN data
CT_eGRN_DAR_Prune_Metadata.l <- readRDS(paste0(RES_DIR, "eRegulons_CT_Prune/CT_eGRN_DAR_Prune_Metadata.rds"))
Peak_propCT.l <- readRDS(paste0(RES_DIR, "eRegulons_CT_Prune/Peak_propPerCT_list.rds"))

# snRNA-seq data
FL.SeuratObj <- readRDS("~/local_data/proj/Dev_Multiome/data/FL_scrna_seurat_20250401.rds")
## Cell Metadata
cell_metadata <- FL.SeuratObj@meta.data

# Set function
COMP_Avg_Prop <- function(GeneSet, FL.SeuratObj, cell_metadata){
  ## Subset Seurat
  sub.SeuratObj <- subset(FL.SeuratObj, features = GeneSet)
  expr.m <- GetAssayData(sub.SeuratObj, layer = "data")
  
  ## Compute AVG Exp
  avg_expr <- expr.m %>%
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
  norm_expr.df <- as.data.frame(zscore_expr)
  norm_expr.df <- norm_expr.df[, CT_ORDER]
  
  ## Compute Exp Proportion
  cell_metadata$cell <- rownames(cell_metadata)
  expr.bool <- expr.m > 0
  Gene_propCT.df <- as_tibble(expr.bool, rownames = "Gene") %>%
    pivot_longer(
      -Gene,
      names_to  = "cell",
      values_to = "expressed"
    ) %>%
    left_join(cell_metadata, by = "cell") %>%
    group_by(Gene, anno_wnn_v51) %>%
    summarise(
      prop_expressed = mean(expressed),
      .groups = "drop"
    )
  
  return(list(sub.SeuratObj = sub.SeuratObj,
              norm_expr.df = norm_expr.df,
              Gene_propCT.df = Gene_propCT.df))
}

##########################
# Which networks to draw #
# Transient
# Cell cycle
# Differential Targets
# Differential Regulators
##########################

################################
################################
# RUNX1 Comparison in HSC & MK #
################################
################################

# Setup
TARGET_CT <- c("HSC", "MK")
TARGET_TF <- "RUNX1"

############################
# Get Differential Targets #
############################

# Target genes
DT.df <- list()
for (ct in TARGET_CT) {
  tmp <- CT_eGRN_DAR_Prune_Metadata.l[[ct]][[1]] %>% 
    filter(TF == TARGET_TF) %>%
    pull(Gene) %>%
    unique() %>%
    sort()
  DT.df[[ct]] <- tmp
}
DT.df[["Share"]] <- intersect(DT.df[[1]], DT.df[[2]])
for (ct in TARGET_CT) {
  DT.df[[ct]] <- setdiff(DT.df[[ct]], DT.df[["Share"]])
}

# Target Regions
DTR.df <- list()
for (ct in TARGET_CT) {
  tmp <- CT_eGRN_DAR_Prune_Metadata.l[[ct]][[1]] %>% 
    filter(TF == TARGET_TF) %>%
    pull(Region) %>%
    unique() %>%
    sort()
  DTR.df[[ct]] <- tmp
}
DTR.df[["Share"]] <- intersect(DTR.df[[1]], DTR.df[[2]])
for (ct in TARGET_CT) {
  DTR.df[[ct]] <- setdiff(DTR.df[[ct]], DTR.df[["Share"]])
}

##########################
# Compute Gene attribute #
##########################

# Get DEGs for CL
DEG.df <- read.csv("/work/DevM_analysis/01.annotation/10.integration_joint_clean/data/FL_wnn_markerGenes.v00.csv")
DEG_Network.df <- DEG.df %>% 
  dplyr::filter(group %in% TARGET_CT) %>%
  dplyr::filter(logfoldchanges > 0.3) %>%
  dplyr::filter(pvals_adj <= 0.05)

# Get Gene Prop
GeneSet.v <- do.call(c, DT.df) %>% unname() %>% as.vector()
COMP_Avg_Prop.o <- COMP_Avg_Prop(GeneSet = GeneSet.v, 
                                 FL.SeuratObj = FL.SeuratObj,
                                 cell_metadata = cell_metadata)
Gene_propCT.df <- COMP_Avg_Prop.o$Gene_propCT.df %>%
  filter(anno_wnn_v51 %in% TARGET_CT)
norm_expr.df <- COMP_Avg_Prop.o$norm_expr.df[, TARGET_CT] %>% 
  mutate(Gene = rownames(.))

#########################################
# Extract graphs with attribute to plot #
#########################################

# Setup
ct <- TARGET_CT[2]
other_ct <- setdiff(TARGET_CT, ct)
CT_eGRN.df <- CT_eGRN_DAR_Prune_Metadata.l[[ct]][[1]] %>%
  filter(TF == TARGET_TF)
other_CT_eGRN.df <- CT_eGRN_DAR_Prune_Metadata.l[[other_ct]][[1]] %>%
  filter(TF == TARGET_TF) %>% 
  filter(Region %in% DTR.df[[other_ct]])
CT_eGRN.df$KEEP <- paste0("keep_in_", ct)
other_CT_eGRN.df$KEEP <- paste0("depleted_in_", ct)

# Create edge
CT_edges <- bind_rows(
  CT_eGRN.df %>% transmute(from = TF, to = Region, 
                           weight = importance_TF2G,
                           class = KEEP),
  CT_eGRN.df %>% transmute(from = Region, to = Gene, 
                           weight = importance_R2G,
                           class = KEEP)
) %>% distinct()
other_CT_edges <- bind_rows(
  other_CT_eGRN.df %>% transmute(from = TF, to = Region, 
                                 weight = importance_TF2G,
                                 class = KEEP),
  other_CT_eGRN.df %>% transmute(from = Region, to = Gene, 
                                 weight = importance_R2G,
                                 class = KEEP)
) %>% distinct()
## Merge edges
edges_fusion <- rbind(CT_edges, other_CT_edges) %>% distinct()
write_tsv(edges_fusion, paste0("./results/CT_eGRN_DAR_Networks/RUNX1_HSC_MK/edges_fusion_RUNX1_in_", ct, ".tsv"))

# Annotate nodes
DEG_ct.df <- DEG_Network.df %>%
  dplyr::filter(group == ct)
nodes_fusion <- tibble(id = unique(c(edges_fusion$from, edges_fusion$to))) %>%
  mutate(
    nodeType = case_when(
      id %in% TARGET_TF ~ "TF",
      id %in% edges_fusion$from ~ "Region",
      id %in% edges_fusion$to ~ "Gene",
      TRUE ~ "Other"
    ),
    nodeClass = case_when(
      id %in% c(TARGET_TF, DTR.df$Share, DT.df$Share) ~ "share",
      id %in% c(DTR.df[[ct]], DT.df[[ct]]) ~ ct,
      id %in% c(DTR.df[[other_ct]], DT.df[[other_ct]]) ~ other_ct,
      TRUE  ~ "Other"
    ),
    nodeDEG = case_when(
      id %in% DEG_ct.df$names ~ "DEG",
      TRUE  ~ "Non-DEG"
    )
  )

# Get Gene Exp & ExpProp
expr_wide <- Gene_propCT.df %>%
  dplyr::rename(cellType = anno_wnn_v51) %>%
  filter(cellType %in% TARGET_CT) %>%
  group_by(Gene, cellType) %>%
  pivot_wider(
    names_from  = cellType,
    values_from = prop_expressed,
    names_prefix = "prop_"
  )
expr_wide <- expr_wide %>% 
  left_join(norm_expr.df, by = "Gene")
nodes_fusion <- nodes_fusion %>%
  left_join(expr_wide, by = c("id" = "Gene"))

# Get Region Openness Prop
region_prop <- Peak_propCT.l[[ct]]
region_prop_df <- enframe(region_prop, name = "id", value = "region_prop") %>%
  mutate(region_prop = as.numeric(region_prop))
nodes_fusion <- nodes_fusion %>%
  left_join(region_prop_df, by = "id")
nodes_fusion$prop_HSC[is.na(nodes_fusion$prop_HSC)] <- 0
nodes_fusion$prop_MK[is.na(nodes_fusion$prop_MK)] <- 0
nodes_fusion$region_prop[is.na(nodes_fusion$region_prop)] <- 0
nodes_fusion$prop_HSC <- nodes_fusion$prop_HSC/max(nodes_fusion$prop_HSC) + 
  nodes_fusion$region_prop/max(nodes_fusion$region_prop)
nodes_fusion$prop_MK <- nodes_fusion$prop_MK/max(nodes_fusion$prop_MK) + 
  nodes_fusion$region_prop/max(nodes_fusion$region_prop)
write_tsv(nodes_fusion, paste0("./results/CT_eGRN_DAR_Networks/RUNX1_HSC_MK/nodes_fusion_RUNX1_in_", ct, ".tsv"))

###############
# GO analysis #
###############

library(clusterProfiler)
library(org.Hs.eg.db)
TG.l <- list()
TG.l[["HSC"]] <- DT.df$HSC
TG.l[["MK"]] <- DT.df$MK
## Covert SYMBOL to ENTREZID
TG_ENTREZID.l <- list()
GOres_list <- list()
for (ct in TARGET_CT) {
  df <- bitr(TG.l[[ct]],
             fromType = "SYMBOL",
             toType = "ENTREZID",
             OrgDb = org.Hs.eg.db)
  TG_ENTREZID.l[[ct]] <- unique(df$ENTREZID)
  ## GO Biological Process enrichment
  eGO.o <- enrichGO(gene = TG_ENTREZID.l[[ct]],
                    OrgDb = org.Hs.eg.db,
                    keyType = "ENTREZID",
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.1,
                    qvalueCutoff  = 0.1, 
                    readable = TRUE)
  ## Simplify GO terms
  simp.o <- simplify(eGO.o,
                     cutoff  = 0.7,
                     by = "p.adjust",
                     select_fun = min)
  GOres_list[[ct]] <- as.data.frame(simp.o@result)
  
  p <- dotplot(simp.o, showCategory = 20) +
    ggtitle(paste0("GO BP of ", ct)) +
    theme_minimal() +
    theme(
      plot.title = element_text(face="bold", hjust = 0.5)
    )
  ggsave(filename = paste0(FIG_DIR, "05.CTeGRN_DAR_Network/RUNX1_HSC_MK/GO_BP_", TARGET_TF, "_", ct, ".pdf"),
         plot = p, width = 6, height = 10)
}
### Write GO terms
library(writexl)
write_xlsx(
  GOres_list,
  path = "./results/CT_eGRN_DAR_Networks/RUNX1_HSC_MK/RUNX1_HSC_MK_GO_list_readable.xlsx"
)

##################
# Compute AUCell #
##################

# Load packages & data
library(AUCell)
library(GSEABase)
## Load snRNA-seq data with WNN UMAP
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
## Load gene rank per cell
cells_rankings <- readRDS("./results/eRegulons_AUCell/GeneRank_for_AUCell.rds")
cells_rankings2 <- cells_rankings[, colnames(FL.SeuratObj)]

#  Geneset setup
geneSets.l <- list()
for (ct in TARGET_CT) {
  geneSets.l[[ct]] <- GeneSet(TG.l[[ct]], setName=as.character(paste0(TARGET_TF, "_", ct)))
}

# Run AUCell
geneSets <- GeneSetCollection(geneSets.l)
names(geneSets)
## Compute gene based AUCell score
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
saveRDS(cells_AUC, paste0("./results/eRegulons_CT_Prune/CT_eGRN_DAR_Prune_AUCell/CT_eGRN_DAR_Prune_", TARGET_TF, "_AUCell.rds"))

# Add new AUCell data to Seurat Object
cells_AUC.m <- assay(cells_AUC)
rownames(cells_AUC.m) <- paste0(TARGET_TF, "_", TARGET_CT)
cells_AUC.m <- cells_AUC.m[, colnames(FL.SeuratObj)]
cells_AUC.df <- as.data.frame(t(cells_AUC.m))
FL.SeuratObj <- AddMetaData(
  object = FL.SeuratObj,
  metadata = cells_AUC.df
)
## Featureplot
idx <- paste0(TARGET_TF, "_", TARGET_CT[2])
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
# p
save_plot(filename = paste0("./plots/05.CTeGRN_DAR_Network/RUNX1_HSC_MK/", idx, "_UMAP.pdf"),
          p, base_height = 10, base_width = 10)





