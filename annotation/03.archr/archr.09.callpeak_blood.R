#!/usr/bin/env Rscript

# run
# nohup Rscript code/archr.10.callpeak_blood.R > logs/archr.10.callpeak_blood.log &

DOCNAME <- "archr_48FL.blood"

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
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "celltype", force = TRUE)

# call peak
proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = "celltype",
  pathToMacs2 = "/work/home/bin/macs2",
  force = TRUE
)

# Add Peak Matrix
proj <- addPeakMatrix(proj, force = TRUE)
getAvailableMatrices(proj)

# save
saveArchRProject(ArchRProj = proj, load = FALSE)

Sys.Date()
sessionInfo()


