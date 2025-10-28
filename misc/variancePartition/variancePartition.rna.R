#!/usr/bin/env Rscript

# cd ~/project/20231127_DevM/devm_r432
# nohup Rscript code/variancePartition.rna.R > logs/variancePartition.rna.log &

DOCNAME <- "variancePartition.rna"
dir.create(here::here("output", DOCNAME), showWarnings = FALSE)

library(tidyverse)
# srat
library(Seurat)
library(Matrix)
# model
library(variancePartition)

# load
srat <- readRDS("/work/DevM_analysis/02.abundance/Milo_FL_PCW250401/data/FL_rna_clustered.rds")
srat

# 1% cell
flag_rna <- rowSums(srat@assays$RNA$counts > 0)/ncol(srat) > 0.01
srat <- srat[flag_rna,]
srat

# modify meta
srat@meta.data <- srat@meta.data |>
  rownames_to_column("barcode") |>
  mutate(PCWsca = as.numeric(scale(as.integer(as.character(PCW))))) |>
  mutate(celltype = anno_wnn_v51,
         Sex = factor(Sex, levels = c("M", "F")),
         donorID = factor(donorID),
         Batch = factor(Batch)) |>
  droplevels() |>
  column_to_rownames("barcode")
str(srat@meta.data)

# get mdata for sample
mdata <- srat@meta.data |>
  rownames_to_column("barcode") |>
  dplyr::select(sampleID, libraryID, donorID, Sex, PCWsca, Batch) |>
  distinct() |>
  droplevels()
str(mdata)
saveRDS(mdata, here::here("output", DOCNAME, "mdata_sample.rds"))

# pseudo-bulk by sample
srat_pb <- AggregateExpression(srat, group.by = "sampleID", return.seurat = TRUE)
mdata_pb <- srat_pb@meta.data |>
  left_join(mdata, by = "sampleID") |>
  mutate(
    sampleID = factor(sampleID)
  ) |>
  droplevels() |>
  column_to_rownames("orig.ident")
str(mdata_pb)

# varPar with PB
mat <- srat_pb@assays$RNA$data
form <- ~ PCWsca + (1 | libraryID) + (1 | donorID) + (1 | Sex) + (1 | Batch)
# Fit model and extract results
varPart <- fitExtractVarPartModel(mat, form, mdata_pb)
# sort variables (i.e. columns) by median fraction of variance explained
vp <- sortCols(varPart)
saveRDS(varPart, file = here::here("output", DOCNAME, "PBvarPart.rds"))
saveRDS(vp, file = here::here("output", DOCNAME, "PBvp.rds"))


sessionInfo()
Sys.Date()
