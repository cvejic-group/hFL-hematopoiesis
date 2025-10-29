#!/usr/bin/env Rscript

.libPaths('/work/home/project/20231127_DevM/devm_rproj/renv/library/R-4.3/x86_64-pc-linux-gnu')

library(tidyverse)

# prep
cellmeta <- read.csv("/work/DevM_analysis/01.annotation/10.integration_joint_clean/data/FL_wnn_cellmeta.v01.csv") %>%
  dplyr::select(barcode = X, libraryID, celltype = anno_wnn_v51) %>%
  mutate(barcode = str_split(barcode, "_", simplify=T)[,2])


# out
out_dir <- "/work/home/project/20231127_DevM/sinto/bam_per_lib"
for (lib in unique(cellmeta$libraryID)) {
  dir.create(file.path(out_dir, lib), showWarnings = FALSE)
  cellmeta %>%
    filter(libraryID == lib) %>%
    dplyr::select(barcode, celltype) %>%
    write_tsv(file.path(out_dir, lib, "cells.tsv"), col_names = FALSE)
}
