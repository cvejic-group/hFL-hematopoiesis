#!/usr/bin/env Rscript

# usage
# Rscript run_scDblFinder.R sample_name

# traceback
#options(error = function() traceback(2))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(scDblFinder))

# args
args <- commandArgs(trailingOnly = TRUE)

# test args
if (length(args) == 0) {
  stop("At least one argument must be supplied (sample name).", call. = FALSE)
} else {
  sam <- args[1]
  cr <- args[2]
}

message(sam)

# input
data_dir <- "~/data/cellranger-arc"

# output dir
work_dir = "~/project/20231127_DevM/scDblFinder"

# build seurat obj (cellranger-arc output)
counts <- Read10X_h5(file.path(data_dir, cr, "filtered_feature_bc_matrix.h5"))
srat <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# rm low UMI cells to avoid error of `Size factors should be positive`
# https://github.com/plger/scDblFinder/issues/96
srat <- subset(srat, subset = nCount_RNA > 100)

# scDblFinder, rm doublet
# seurat v4
set.seed(123)
srat <- as.Seurat(scDblFinder(as.SingleCellExperiment(srat)))

# save
srat@meta.data %>%
  rownames_to_column("barcode") %>%
  dplyr::select(barcode, scDblFinder.class, scDblFinder.score, scDblFinder.weighted, scDblFinder.cxds_score) %>%
  write_csv(file = file.path(work_dir, paste0(sam, ".csv")))

# session info
devtools::session_info()
