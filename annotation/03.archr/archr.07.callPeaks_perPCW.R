#!/usr/bin/env Rscript

# run
# nohup Rscript analysis/archr.08.callPeaks_perPCW.R > analysis/archr.08.callPeaks_perPCW.log &

DOCNAME <- "archr_48FL"
dir.create(here::here("output", DOCNAME, "MarkerFeatures"), showWarnings = FALSE)

library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyverse)
library(ArchR)
#reticulate::use_python("/work/home/software/anaconda3/bin/python")

# set up
set.seed(1)
# default number of threads
addArchRThreads(threads = 16)
# add a reference genome annotation for ArchR
addArchRGenome("hg38")
# chr prefix to exclude KI/GL scaffolds
addArchRChrPrefix(chrPrefix = TRUE)
# cor cutoff
corCutOff = 0.5

# load
proj <- loadArchRProject(path = here::here('output', DOCNAME))
proj

# addGroupCoverages
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "PCW_celltype", force = TRUE)

# call peak
proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = "PCW_celltype",
  pathToMacs2 = "/work/home/bin/macs2",
  force = TRUE
)

# save
saveArchRProject(ArchRProj = proj, load = FALSE)

