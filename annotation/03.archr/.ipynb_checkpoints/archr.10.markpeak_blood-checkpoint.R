#!/usr/bin/env Rscript

# run
# nohup Rscript code/archr.11.markpeak_blood.R > logs/archr.11.markpeak_blood.log &

DOCNAME <- "archr_48FL.blood"
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

# marker peak
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "celltype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
saveRDS(markersPeaks, file = here::here('output', DOCNAME, "MarkerFeatures", "markersPeaks.rds"))


Sys.Date()
sessionInfo()


