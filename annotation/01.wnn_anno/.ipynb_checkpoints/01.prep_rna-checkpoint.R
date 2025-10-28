#!/usr/bin/env Rscript

# run
# nohup Rscript analysis/intg_rna_48FL.00.harmony.R > analysis/intg_rna_48FL.00.harmony.log &

DOCNAME <- "intg_rna_48FL.00"
dir.create(here::here("output", DOCNAME), showWarnings = FALSE)

library(tidyverse)

# srat
library(Seurat)
library(SeuratWrappers)
library(harmony)
reticulate::use_condaenv("/work/home/software/anaconda3/envs/scanpy")

source(here::here("code/preprocess.R"))
source(here::here("code/color.R"))
source(here::here("code/plot.R"))

# load library info
sam_info <- read_tsv("/work/DevM_analysis/utils/sam_info.DevM.12.08.24.txt") |>
  dplyr::select(libraryID, donorID, Tissue, PCW, Sex)

# load donor info
donor_info <- sam_info |>
  filter(Sex %in% c("M", "F")) |>
  dplyr::select(-libraryID) |>
  distinct()

# load clean srat obj
data_dir <- "/work/DevM_analysis/01.annotation/02.cleandata/data/seu_obj"
srat_lst <- lapply(unique(sam_info$libraryID), function(i) {
  cat(i)
  fname <- file.path(data_dir, paste0(i, "_rna_seu.rds"))
  srat <- readRDS(fname)
  DefaultAssay(srat) <- "RNA"
  # assay 5 to 3
  srat[["RNA"]] <- as(object = srat[["RNA"]], Class = "Assay")
  # doublet score
  df_doublet <- read_csv(file.path("/work/DevM_analysis/01.annotation/01.doublet/data", paste0(i, ".csv"))) |>
    dplyr::select(barcode, scDblFinder.score)
  # neat metadata
  srat@meta.data <- srat@meta.data |>
    rownames_to_column("barcode") |>
    left_join(df_doublet, by = "barcode") |>
    dplyr::select(
      barcode, nCount_RNA, nFeature_RNA, percent.mt, percent.rb,
      scDblFinder.score,
      libraryID, donorID, sampleID
    ) |>
    column_to_rownames("barcode")
  # rename cells
  srat <- RenameCells(srat, new.names = paste0(i, "_", colnames(srat)))
  return(srat)
})
names(srat_lst) <- unique(sam_info$libraryID)
lapply(srat_lst, dim)

# HVGs in each lib
srat_lst <- lapply(srat_lst, function(i) {
  i |>
    NormalizeData(verbose = FALSE) |>
    FindVariableFeatures(nfeatures = 2000, verbose = FALSE)
})

# merge
# srat <- purrr::reduce(srat_lst, merge)
srat <- merge(srat_lst[[1]], srat_lst[2:length(srat_lst)])
srat
table(srat$libraryID, useNA = "ifany")
table(srat$donorID, useNA = "ifany")

# factor level
srat$libraryID <- factor(srat$libraryID, levels = unique(sam_info$libraryID))
srat$donorID <- factor(srat$donorID, levels = donor_info$donorID)

# donor info (Tissue, PCW, Sex)
srat@meta.data <- srat@meta.data |>
  rownames_to_column("barcode") |>
  left_join(donor_info, by = "donorID") |>
  column_to_rownames("barcode")
table(srat$PCW, useNA = "ifany")
table(srat$Sex, useNA = "ifany")

# HVGs
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
  NormalizeData() |>
  ScaleData(features = VariableFeatures(srat), vars.to.regress = c("nCount_RNA", "nFeature_RNA")) |>
  RunPCA(features = VariableFeatures(srat), npcs = 100)

# cell cycle
srat <- CellCycleScoring(srat,
                         s.features = cc.genes.updated.2019$s.genes,
                         g2m.features = cc.genes.updated.2019$g2m.genes,
                         set.ident = FALSE)

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
saveRDS(srat, file = here::here("output", DOCNAME, "FL_rna_clustered.v00.rds"))
srat@meta.data |>
  rownames_to_column("barcode") |>
  dplyr::select(-starts_with("RNA_snn_res."), -seurat_clusters) |>
  write_csv(file = here::here("output", DOCNAME, "FL_rna_cellmeta.v00.csv"))

# DGE
Idents(srat) <- "RNA_snn_res.3"
marker_df <- FindAllMarkers(srat, only.pos = TRUE)
marker_df |>
  write_csv(here::here("output", DOCNAME, "FL_rna_markerGenes.v00.csv"))
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
  write_csv(here::here("output", DOCNAME, "FL_rna_markerGenes_top100.v00.csv"))

# clean barcodes per lib for ArchR
dir.create(here::here("output", DOCNAME, "clean_barcodes_per_lib"), showWarnings = FALSE)
df <- read.table("/work/DevM_analysis/01.annotation/02.cleandata/data/cell-info_filtering.txt") |>
  filter(HighQualityCell == 1)
barcode_lst <- split(df$gex_barcode, df$libraryID)
tmp <- lapply(1:length(barcode_lst), function(i){
  s = names(barcode_lst)[i]
  print(s)
  data.frame(barcode=barcode_lst[[i]], cell_id=barcode_lst[[i]]) |>
    write_csv(here::here("output", DOCNAME, "clean_barcodes_per_lib", paste0(s, ".csv")))
})

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
                      outFile=here::here("output", DOCNAME, "FL_rna_clustered.v00.h5ad"))

# session info
Sys.Date()
sessionInfo()

