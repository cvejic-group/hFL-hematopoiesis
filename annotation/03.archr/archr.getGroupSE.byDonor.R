#!/usr/bin/env Rscript

# run
# nohup Rscript code/archr.getGroupSE.byDonor.R > logs/archr.getGroupSE.byDonor.log &

DOCNAME <- "archr_48FL"
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
addArchRThreads(threads = 8)
# add a reference genome annotation for ArchR
addArchRGenome("hg38")
# chr prefix to exclude KI/GL scaffolds
addArchRChrPrefix(chrPrefix = TRUE)
# cor cutoff
corCutOff = 0.5

# load
proj <- loadArchRProject(path = here::here('output', DOCNAME))
proj
getAvailableMatrices(proj)

# group by sample
se <- getGroupSE(ArchRProj = proj, useMatrix = "PeakMatrix", groupBy = "donorID", divideN = FALSE, scaleTo = NULL)
saveRDS(se, file = file.path(se_dir, "PeakMatrix.PB_donorID.rds"))

sessionInfo()
Sys.Date()
