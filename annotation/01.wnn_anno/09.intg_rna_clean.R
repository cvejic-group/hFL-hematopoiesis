#!/usr/bin/env Rscript

# run on Rstudio server 432
# nohup Rscript 01.intg_rna_clean.R > 01.intg_rna_clean.log &


.libPaths('/work/home/project/20231127_DevM/devm_rproj/renv/library/R-4.3/x86_64-pc-linux-gnu')

library(tidyverse)

# srat
library(Seurat)
library(SeuratWrappers)
library(harmony)
reticulate::use_condaenv("/work/home/software/anaconda3/envs/scanpy")
source("/work/DevM_analysis/utils/utils.R")

# anno
new_anno = "anno_wnn_v5"
df_anno = read_csv("/work/DevM_analysis/01.annotation/09.annotation_joint/data/FL_wnn_cellmeta.v05.csv") |>
  dplyr::select(barcode = `...1`, !!sym(new_anno), libraryID)

# subset
work_dir = "/work/DevM_analysis/01.annotation/10.integration_joint_clean"
celltype_to_sub = "LowQ"
data_dir = file.path(work_dir, "data")
srat <- readRDS("/work/DevM_analysis/01.annotation/05.integration_rna/data/FL_rna_clustered.v00.rds")
srat@meta.data <- srat@meta.data |>
  rownames_to_column("barcode") |>
  left_join(df_anno, by = "barcode") |>
  column_to_rownames("barcode")
srat <- subset(srat, subset = anno_wnn_v5 == celltype_to_sub, invert = TRUE)
srat <- DietSeurat(srat, layers = c("counts", "data"))
srat

# levels
srat$anno_wnn_v5 <- factor(srat$anno_wnn_v5, levels = c("HSC", "35-MDP?", "GP", "Granulocyte",
                  "49-MEMP", "MEMP", "MEP", "MEMP-Mast-Ery", "MEMP-Ery", "Early-Ery", "Late-Ery",
                  "MEMP-MK", "MK", "49:2-MastP", "MastP", "Mast",
                  "Monocyte", "Kupffer", "cDC1", "cDC2", "pDC", "ASDC",
                  "LMPP", "LP", "Cycling_LP", "PreProB", "ProB-1", "ProB-2", "Large-PreB", "Small-PreB", "IM-B",
                  "NK", "ILCP", "T",
                  "Hepatocyte", "Endothelia"))
table(srat$anno_wnn_v5, useNA = "ifany")


# HVGs in each lib
srat_lst <- SplitObject(srat, split.by = "libraryID")
srat_lst <- lapply(srat_lst, function(i) {
  i |>
    FindVariableFeatures(nfeatures = 2000, verbose = FALSE)
})
HVGs <- SelectIntegrationFeatures(
  object.list = srat_lst,
  nfeatures = 3000
)
length(HVGs)
length(HVGs[grep("MT-", HVGs)])
length(HVGs[grep("^RP[SL]", HVGs)])
# remove MT from HVGs
VariableFeatures(srat) <- grep(pattern = "MT-", x = HVGs, value = TRUE, invert = TRUE)
length(VariableFeatures(srat))

# PCA
srat <- srat |>
  ScaleData(features = VariableFeatures(srat), vars.to.regress = c("nCount_RNA", "nFeature_RNA")) |>
  RunPCA(features = VariableFeatures(srat), npcs = 100)

# PC number
ncomps <- ncomps_signi_pegasus(srat)
srat@misc$ncomps <- ncomps
ncomps

# Harmony theta=2 for library and donor
set.seed(1)
srat <- RunHarmony(srat,
                   group.by.vars = c("libraryID", "donorID"),
                   theta = c(2, 2),
                   lambda = c(1, 1),
                   dims.use = 1:ncomps,
                   plot_convergence = TRUE
)
srat <- RunUMAP(srat, n.neighbors = 20, reduction = "harmony", dims = 1:ncomps)
srat <- FindNeighbors(srat, reduction = "harmony", dims = 1:ncomps)
srat <- FindClusters(srat,
                     resolution = c(0.5, 1, 2, 3, 4, 5, 6),
                     algorithm = 4, method = "igraph"
)

# save
srat
saveRDS(srat, file = file.path(data_dir, "FL_rna_clustered.v00.rds"))
srat@meta.data |>
  rownames_to_column("barcode") |>
  write_csv(file = file.path(data_dir, "FL_rna_cellmeta.v00.csv"))

# to h5ad
new_obj <- srat
# rm useless columns
new_obj@meta.data <- new_obj@meta.data |>
  rownames_to_column("barcode") |>
  dplyr::select(-starts_with("RNA_snn_res."), -seurat_clusters) |>
  column_to_rownames("barcode")
# assay5 to assay3 (sceasy works only for old Seurat obj)
#new_obj[["RNA"]] <- as(object = new_obj[["RNA"]], Class = "Assay")
new_obj
# HVGs compatible with scanpy
new_obj[['RNA']]@meta.features = data.frame(row.names = rownames(new_obj),
                                            highly_variable = ifelse(rownames(new_obj) %in% VariableFeatures(new_obj),
                                                                     TRUE, FALSE))
sceasy::convertFormat(new_obj, from="seurat", to="anndata",
                      main_layer = "counts",
                      drop_single_values=FALSE,
                      outFile=file.path(data_dir, "FL_rna_clustered.v00.h5ad"))

# DGE for RNA anno (the resolution is based on previous observation)
Idents(srat) <- "RNA_snn_res.3"
marker_df <- FindAllMarkers(srat, only.pos = TRUE)
marker_df |>
  write_csv(file.path(data_dir, "FL_rna_markerGenes.v00.csv"))
marker_df_top <- marker_df |>
  group_by(cluster) |>
  dplyr::filter(avg_log2FC > 1)  |>
  slice_head(n = 100) |>
  ungroup()
l <- split(marker_df_top$gene, marker_df_top$cluster)
marker_df_top <- stringi::stri_list2matrix(l, byrow=FALSE) |>
  as.data.frame()
colnames(marker_df_top) <- paste0("cluster_", 1:length(l))
marker_df_top |>
  write_csv(file.path(data_dir, "FL_rna_markerGenes_top100.v00.csv"))

# clean barcodes per lib for ArchR
dir.create(file.path(data_dir, "clean_barcodes_per_lib"), showWarnings = FALSE)
df <- df_anno |>
  filter(anno_wnn_v5 != "LowQ")
barcode_lst <- split(df$barcode, df$libraryID)
tmp <- lapply(1:length(barcode_lst), function(i){
  s = names(barcode_lst)[i]
  print(s)
  data.frame(barcode=barcode_lst[[i]], cell_id=barcode_lst[[i]]) |>
    mutate(barcode=str_split(barcode, '_', simplify = TRUE)[,2],
          cell_id=str_split(cell_id, '_', simplify = TRUE)[,2]
          ) |>
    write_csv(
        file.path(data_dir, "clean_barcodes_per_lib", paste0(s, ".csv")))
})


# session info
Sys.Date()
sessionInfo()


