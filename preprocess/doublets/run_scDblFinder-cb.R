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
}

# sam = "M16216-A"
message(sam)

# input
data_dir <- "~/project/20231127_DevM/cellbender"

# output dir
work_dir = "~/project/20231127_DevM/scDblFinder"
out_dir <- file.path(work_dir, sam)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# build seurat obj (cellbender output)
counts <- Read10X_h5(file.path(data_dir, sam, "output_filtered_srat.h5"))
srat <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)
# scDblFinder, rm doublet
# seurat v4
srat <- as.Seurat(scDblFinder(as.SingleCellExperiment(srat)))

# save
srat@meta.data %>%
  rownames_to_column("barcode") %>%
  dplyr::select(barcode, scDblFinder.class, scDblFinder.score, scDblFinder.weighted, scDblFinder.cxds_score) %>%
  write_csv(file = file.path(out_dir, "scDblFinder_output.csv"))


