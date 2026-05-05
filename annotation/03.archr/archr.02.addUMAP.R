#!/usr/bin/env Rscript

# run
# nohup Rscript analysis/archr.02.addUMAP.R > analysis/archr.02.addUMAP.log &

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

# LSI - ATAC
proj <- addIterativeLSI(
  ArchRProj = proj,
  useMatrix = "TileMatrix",
  depthCol = "nFrags", # default
  name = "LSI_ATAC",
  iterations = 2, # default
  #clusterParams = list(resolution = c(1), maxClusters = 20),
  varFeatures = 25000, # large value possible cause p[length(p)] cannot exceed 2^31-1
  dimsToUse = 1:30, # base on the benchmark, 10-50 is fine
  corCutOff = corCutOff,
  force = TRUE
)

# harmony based on ATAC
proj <- addHarmony(
  ArchRProj = proj,
  reducedDims = "LSI_ATAC",
  name = "Harmony",
  groupBy = c("libraryID", "donorID"),
  corCutOff = corCutOff,
  force = TRUE
)

# umap based on ATAC
proj <- addUMAP(ArchRProj = proj, reducedDims = "Harmony",
                name = "UMAP",
                corCutOff = corCutOff,
                force = TRUE
)

# impute gene score based on ATAC
proj <- addImputeWeights(proj, reducedDims = "Harmony", corCutOff = corCutOff)

# save
saveArchRProject(ArchRProj = proj, load = FALSE)

Sys.Date()
sessionInfo()


