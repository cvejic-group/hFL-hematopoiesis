#!/usr/bin/env Rscript

# run
# nohup Rscript analysis/archr.03.getMarkerFeatures.R > analysis/archr.03.getMarkerFeatures.log &

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

# marker genes by anno_wnn_v5
markersGS <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  groupBy = "celltype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
saveRDS(markersGS, file = here::here('output', DOCNAME, "MarkerFeatures", "markersGS.rds"))

# addGroupCoverages
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "celltype")

# call peak
proj <- addReproduciblePeakSet(
  ArchRProj = proj,
  groupBy = "celltype",
  pathToMacs2 = "/work/home/bin/macs2"
)

# Add Peak Matrix
proj <- addPeakMatrix(proj)
getAvailableMatrices(proj)

# marker peak
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "celltype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
saveRDS(markersPeaks, file = here::here('output', DOCNAME, "MarkerFeatures", "markersPeaks.rds"))

# add motif
if("Motif" %ni% names(proj@peakAnnotation)){
  proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")
}
# motif enrichment
enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
saveRDS(enrichMotifs, file = here::here('output', DOCNAME, "MarkerFeatures", "enrichMotifs.rds"))

# save
saveArchRProject(ArchRProj = proj, load = FALSE)

Sys.Date()
sessionInfo()

