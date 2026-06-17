#!/usr/bin/env Rscript

# run on Rstudio server 432

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
  dplyr::select(barcode = `...1`, celltype = !!sym(new_anno))

# subset
work_dir = "/work/DevM_analysis/01.annotation/11.subclustering"
celltype_to_sub = "HSC"
data_dir = file.path(work_dir, celltype_to_sub, "data")
srat <- readRDS("/work/DevM_analysis/01.annotation/05.integration_rna/data/FL_rna_clustered.v00.rds")
srat@meta.data <- srat@meta.data |>
  rownames_to_column("barcode") |>
  left_join(df_anno, by = "barcode") |>
  column_to_rownames("barcode")
srat <- subset(srat, subset = celltype == celltype_to_sub)
srat <- DietSeurat(srat, layers = c("counts", "data"))

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

# session info
Sys.Date()
sessionInfo()


