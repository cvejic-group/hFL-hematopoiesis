#!/usr/bin/env Rscript

# nohup Rscript code/07.ROGUE.R > logs/07.ROGUE.log &

DOCNAME <- "07.ROGUE"
dir.create(here::here("output", DOCNAME), showWarnings = FALSE)

library(tidyverse)
library(Seurat)
library(ROGUE)

# load
df <- read.csv("/work/DevM_analysis/01.annotation/11.subclustering/HSC/data/FL_wnn_cellmeta.v00.csv") |>
  dplyr::select(barcode = X, starts_with("leiden_wnn_"))
srat <- readRDS("/work/DevM_analysis/01.annotation/11.subclustering/HSC/data/FL_rna_clustered.v00.rds")
srat@meta.data <- srat@meta.data |>
  rownames_to_column("barcode") |>
  left_join(df, by = "barcode") |>
  mutate(LibraryDonor = paste(libraryID, donorID, sep = "@")) |>
  column_to_rownames("barcode")
expr <- GetAssayData(srat, layer = "counts")

# functions
cal_rogue <- function(sratObj = NULL, cluster = NULL) {
  expr <- GetAssayData(sratObj, layer = "counts")
  ent.res <- SE_fun(expr)
  rogue.value <- CalculateRogue(ent.res, platform = "UMI")
  title <- paste(
    paste0("Cluster ", cluster),
    paste0("(ROGUE = ", rogue.value, ")")
  )
  p <- SEplot(ent.res) +
    labs(title = title)
  return(p)
}
cal_rogue_all <- function(sratObj = NULL, cluster_key = NULL) {
  cluster_num <- max(as.numeric(sratObj@meta.data[[cluster_key]]))
  sratObj@meta.data[[cluster_key]] <- factor(sratObj@meta.data[[cluster_key]],
                                             levels = as.character(0:cluster_num))
  plst <- lapply(levels(sratObj@meta.data[[cluster_key]]), function(x) {
    srat_sub <- subset(sratObj, subset = !!sym(cluster_key) == x)
    cal_rogue(sratObj = srat_sub, cluster = x)
  })
  return(plst)
}

# leiden_wnn_0.1
cluster_key <- "leiden_wnn_0.1"

# by this way we get different ROGUE values from the ones of each sample, so discard this method
#plst <- cal_rogue_all(sratObj = srat, cluster_key = cluster_key)
#saveRDS(plst, file = here::here("output", DOCNAME, paste0(cluster_key, ".plst.rds")))

# low cell number generate low ROGUE values (which cannot be trusted?), so I use min.cell.n = 100 here
rogue.res <- rogue(expr, labels = srat@meta.data[[cluster_key]], samples = srat$sampleID, platform = "UMI", span = .6, min.cell.n = 100)
saveRDS(rogue.res, file = here::here("output", DOCNAME, paste0(cluster_key, ".rogue.rds")))

# leiden_wnn_0.3
cluster_key <- "leiden_wnn_0.3"
rogue.res <- rogue(expr, labels = srat@meta.data[[cluster_key]], samples = srat$sampleID, platform = "UMI", span = .6, min.cell.n = 100)
saveRDS(rogue.res, file = here::here("output", DOCNAME, paste0(cluster_key, ".rogue.rds")))

# leiden_wnn_0.5
cluster_key <- "leiden_wnn_0.5"
rogue.res <- rogue(expr, labels = srat@meta.data[[cluster_key]], samples = srat$sampleID, platform = "UMI", span = .6, min.cell.n = 100)
saveRDS(rogue.res, file = here::here("output", DOCNAME, paste0(cluster_key, ".rogue.rds")))

# leiden_wnn_1
cluster_key <- "leiden_wnn_1"
rogue.res <- rogue(expr, labels = srat@meta.data[[cluster_key]], samples = srat$sampleID, platform = "UMI", span = .6, min.cell.n = 100)
saveRDS(rogue.res, file = here::here("output", DOCNAME, paste0(cluster_key, ".rogue.rds")))

Sys.Date()
sessionInfo()
