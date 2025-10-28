#!/usr/bin/env Rscript

# cd ~/project/20231127_DevM/devm_r432
# nohup Rscript code/variancePartition.atac.R > logs/variancePartition.atac.log &

DOCNAME <- "variancePartition.atac"
dir.create(here::here("output", DOCNAME), showWarnings = FALSE)

library(tidyverse)
# srat
library(SummarizedExperiment)
library(SingleCellExperiment)
library(scater)
library(Seurat)
library(Matrix)
# model
library(variancePartition)

# load
se <- readRDS("~/project/20231127_DevM/archr_rproj/output/archr_48FL.blood/GroupSE/PeakMatrix.PB_sampleID.rds")
se <- as(se, "SingleCellExperiment")
assays(se)$counts <- assays(se)$PeakMatrix
se <- logNormCounts(se)
srat_pb <- as.Seurat(se)
srat_pb <- RenameAssays(srat_pb, originalexp = "RNA")
srat_pb

# 1% cell
peak_fraction <- read_csv("~/project/20231127_DevM/archr_rproj/output/archr_48FL.blood/GroupSE/PeakFractionInCells.tsv")
flag_atac <- peak_fraction$Fraction > 0.01
srat_pb <- srat_pb[flag_atac,]
srat_pb

# get mdata for sample
mdata_pb <- readRDS(here::here("output", "variancePartition.rna", "mdata_sample.rds")) |>
  column_to_rownames("sampleID")
srat_pb <- srat_pb[, rownames(mdata_pb)]
srat_pb

# modify mdata
mdata_pb <- mdata_pb[colnames(srat_pb),]
dim(mdata_pb)

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

