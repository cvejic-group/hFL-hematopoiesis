#!/usr/bin/env Rscript

# run
# nohup Rscript code/05.get_atac_pb_by_donor_promoter.R

DOCNAME <- "archr_48FL.pro_peak_mat"
# dirs
archr_dir <- here::here("output", DOCNAME)
se_dir <- file.path(archr_dir, "GroupSE")
dir.create(se_dir, showWarnings = FALSE, recursive = TRUE)

library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyverse)
library(ArchR)
#reticulate::use_python("/work/home/software/anaconda3/bin/python")

# set up
set.seed(1)
# default number of threads
addArchRThreads(threads = 4)
# add a reference genome annotation for ArchR
addArchRGenome("hg38")
# chr prefix to exclude KI/GL scaffolds
addArchRChrPrefix(chrPrefix = TRUE)
# cor cutoff
corCutOff = 0.5

# load
proj <- loadArchRProject(path = here::here('output', DOCNAME))
proj

# gt_pro
gr_pro <- readRDS("/work/home/project/20231127_DevM/devm_r432/data/cellranger-2020-A-2.0.0.promoter.rds")

# check
getAvailableMatrices(ArchRProj = proj)

# add peak
proj <- addPeakSet(
  ArchRProj = proj,
  peakSet = gr_pro,
  force = TRUE
)

# add peak mat
proj <- addPeakMatrix(proj)

# save
saveArchRProject(ArchRProj = proj, load = FALSE)

# group by donor
se <- getGroupSE(ArchRProj = proj, useMatrix = "PeakMatrix", groupBy = "donorID", divideN = FALSE, scaleTo = NULL)
saveRDS(se, file = file.path(se_dir, "PeakMatrix.PB_donorID.rds"))

Sys.Date()
sessionInfo()

