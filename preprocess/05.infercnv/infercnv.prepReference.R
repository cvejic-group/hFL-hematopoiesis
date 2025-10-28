#!/usr/bin/env Rscript

# usage

library(tidyverse)
library(Seurat)

sam_info <- read_tsv(here::here("data/sam_info.tsv")) |>
  filter(SampleSource == "HDBR", !(PCW %in% c(6, 7))) |> # no pooled lib
  dplyr::select(libraryID, donorID)

seu_dir <- "/work/DevM_analysis/01.annotation/02.cleandata/data/seu_obj"
# merge srat for each donor, then subset 500 cells from each donor
srat_ctrl_lst <- lapply(split(sam_info$libraryID, sam_info$donorID), function(sams) {
  srat_lst <- lapply(sams, function(d) {
    srat <- readRDS(file.path(seu_dir, paste0(d, "_rna_seu.rds")))
    srat[["RNA"]] <- as(object = srat[["RNA"]], Class = "Assay")
    srat <- RenameCells(srat, new.names = paste0(d, "_", colnames(srat)))
    return(srat)
  })
  srat_donor <- purrr::reduce(srat_lst, merge)
  Idents(srat_donor) <- "donorID"
  return(subset(srat_donor, downsample = 500, seed = 1234))
})

# merge these down sampled donors
srat_ctrl <- purrr::reduce(srat_ctrl_lst, merge)
srat_ctrl <- DietSeurat(srat_ctrl)
srat_ctrl

table(srat_ctrl$donorID, useNA = "ifany")
table(srat_ctrl$libraryID, useNA = "ifany")

# out
saveRDS(srat_ctrl, file = here::here("data", "mergedHDBRCtrl.4_infercnv.rds"))

