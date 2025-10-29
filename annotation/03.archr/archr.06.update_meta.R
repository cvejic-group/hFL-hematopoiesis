#!/usr/bin/env Rscript

# run
# nohup Rscript analysis/archr.01.ArchRProject.R > analysis/archr.01.ArchRProject.log &

DOCNAME <- "archr_48FL"
dir.create(here::here("output", DOCNAME), showWarnings = FALSE)

library(tidyverse)
library(ArchR)
library(harmony)
#reticulate::use_condaenv("~/software/anaconda3/envs/scvi-env")

# set up
set.seed(1)
# default number of threads
addArchRThreads(threads = 12)
# add a reference genome annotation for ArchR
addArchRGenome("hg38")
# chr prefix to exclude KI/GL scaffolds
addArchRChrPrefix(chrPrefix = TRUE)
# cor cutoff
corCutOff = 0.5

# load
proj <- loadArchRProject(path = here::here('output', DOCNAME))
proj

# lib and donor info
lib_info <- read_csv("/work/DevM_analysis/utils/sam_info.DevM.22.10.24.csv") |>
  filter(Tissue == "FL", Sex != "Unknown")
df_pcw_batch <- lib_info |>
  dplyr::select(libraryID, PCW, Batch) |>
  distinct()

# update current
cellmeta <- getCellColData(proj) %>%
  as.data.frame() %>%
  rownames_to_column("barcode") %>%
  dplyr::select(everything(), -PCW, anno_wnn_v5 = celltype) %>%
  left_join(df_pcw_batch, by = "libraryID") %>% # PCW and Batch
  mutate(
    libraryID = factor(libraryID, levels = df_pcw_batch$libraryID),
    donorID = factor(donorID, levels = unique(lib_info$donorID)),
    sampleID = factor(sampleID, levels = lib_info$sampleID),
    PCW = factor(PCW, levels = c(5:18))
  ) %>%
  mutate(
    anno_wnn_v51 = case_when(
      anno_wnn_v5 == "35-MDP" ~ "MDP",
      anno_wnn_v5 == "49-MEMP" ~ "MEMP-t",
      anno_wnn_v5 == "49:2-MastP" ~ "MastP-t",
      anno_wnn_v5 == "Cycling_LP" ~ "Cycling-LP",
      TRUE ~ anno_wnn_v5
    )
  )
cellmeta$anno_wnn_v51 <- factor(cellmeta$anno_wnn_v51, levels = c("HSC", "GP", "Granulocyte",
                                                          "MEMP-t", "MEMP", "MEP", "MEMP-Mast-Ery", "MEMP-Ery", "Early-Ery", "Late-Ery",
                                                          "MEMP-MK", "MK", "MastP-t", "MastP", "Mast",
                                                          "MDP", "Monocyte", "Kupffer", "cDC1", "cDC2", "pDC", "ASDC",
                                                          "LMPP", "LP", "Cycling-LP", "PreProB", "ProB-1", "ProB-2", "Large-PreB", "Small-PreB", "IM-B",
                                                          "NK", "ILCP", "T",
                                                          "Hepatocyte", "Endothelia"))

table(cellmeta$anno_wnn_v51, useNA = "ifany")

# update
proj$PCW <- cellmeta$PCW
proj$Batch <- cellmeta$Batch
proj$anno_wnn_v5 <- proj$celltype
proj$anno_wnn_v51 <- cellmeta$anno_wnn_v51
proj$celltype <- proj$anno_wnn_v51 # update the celltype column to the latest

# PCW_celltype
proj$PCW_celltype <- paste0( proj$celltype, "_PCW", proj$PCW)

# save
saveArchRProject(ArchRProj = proj, load = FALSE)

