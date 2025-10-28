#!/usr/bin/env Rscript

# run on Rstudio server 432
# nohup Rscript 01.prep_rna.R > 01.prep_rna.log &

.libPaths('/work/home/project/20231127_DevM/devm_rproj/renv/library/R-4.3/x86_64-pc-linux-gnu')

library(tidyverse)

# srat
library(Seurat)
library(SeuratWrappers)
library(harmony)
reticulate::use_condaenv("/work/home/software/anaconda3/envs/scanpy")
source("/work/DevM_analysis/utils/utils.R")

# load
work_dir = "/work/DevM_analysis/02.abundance/Milo_FL_PCW250401"
data_dir = file.path(work_dir, "data")
srat <- readRDS("/work/DevM_analysis/01.annotation/11.subclustering/blood/data/FL_rna_clustered.v00.rds")

# subset - PCW
srat <- subset(srat, subset = PCW != "Mixed")
srat <- DietSeurat(srat, layers = c("counts", "data"))
srat@meta.data <- droplevels(srat@meta.data)
srat

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

# save
srat
saveRDS(srat, file = file.path(data_dir, "FL_rna_clustered.rds"))

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
                      outFile=file.path(data_dir, "FL_rna_clustered.h5ad"))

# session info
Sys.Date()
sessionInfo()


