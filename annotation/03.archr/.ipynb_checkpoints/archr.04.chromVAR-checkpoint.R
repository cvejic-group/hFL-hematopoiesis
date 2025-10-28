#!/usr/bin/env Rscript

# run
# nohup Rscript analysis/archr.04.chromVAR.R > analysis/archr.04.chromVAR.log &

DOCNAME <- "archr_48FL"
dir.create(here::here("output", DOCNAME), showWarnings = FALSE)

library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyverse)
library(ArchR)
#reticulate::use_python("/work/home/software/anaconda3/bin/python")

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
getAvailableMatrices(proj)

# add motif
if("Motif" %ni% names(proj@peakAnnotation)){
  proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")
}

# chromVAR, slow
proj <- addBgdPeaks(proj, force = TRUE)
proj <- addDeviationsMatrix(
  ArchRProj = proj,
  peakAnnotation = "Motif",
  force = TRUE
)
getAvailableMatrices(proj)

# save
saveArchRProject(ArchRProj = proj, outputDirectory = here::here("output", DOCNAME), load = FALSE)

# diff
diffVariability <- getMarkerFeatures(
  ArchRProj = proj,
  testMethod = "wilcoxon",
  binarize = FALSE,
  useMatrix = "MotifMatrix",
  groupBy = "celltype",
  useSeqnames="z"
)
saveRDS(diffVariability, file = here::here('output', DOCNAME, "MarkerFeatures", "diffVariability.rds"))

diffDeviation <- getMarkerFeatures(
  ArchRProj = proj,
  testMethod = "wilcoxon",
  binarize = FALSE,
  useMatrix = "MotifMatrix",
  groupBy = "celltype",
  useSeqnames="deviations"
)
saveRDS(diffDeviation, file = here::here('output', DOCNAME, "MarkerFeatures", "diffDeviation.rds"))

# save
saveArchRProject(ArchRProj = proj, outputDirectory = here::here("output", DOCNAME), load = FALSE)

# save to another one (old as backup)
# the ArchR save sucks (random unknown error happens when reading)
#saveArchRProject(ArchRProj = proj, outputDirectory = here::here("output", "archr_48FL.backup"), load = FALSE)
R.utils::copyDirectory(here::here("output", "archr_48FL"), here::here("output", "archr_48FL.backup"))


Sys.Date()
sessionInfo()

