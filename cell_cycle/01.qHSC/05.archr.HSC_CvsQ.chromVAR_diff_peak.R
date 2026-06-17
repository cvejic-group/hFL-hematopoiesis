#!/usr/bin/env Rscript

# run
# nohup Rscript code/archr.HSC_CvsQ.chromVAR_diff_peak.R > logs/archr.HSC_CvsQ.chromVAR_diff_peak.log &

DOCNAME <- "archr_48FL"

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

# re-add peak mat
peaks <- readRDS("/work/home/project/20231127_DevM/gwas/sLDSC/atac_peaks/archr_peakSet.rds")
proj <- addPeakSet(ArchRProj = proj, peakSet = peaks, force = TRUE)
proj <- addPeakMatrix(proj, force = TRUE)

# HSC_CvsQ diff. peaks
HSC_CvsQ <- c(
  G0up = "/work/home/project/20231127_DevM/HSC.CvsQ.peaks/G0_up.bed",
  G0dn = "/work/home/project/20231127_DevM/HSC.CvsQ.peaks/G0_dn.bed"
)

if("HSC_CvsQ" %ni% names(proj@peakAnnotation)){
  proj <- addPeakAnnotations(ArchRProj = proj, regions = HSC_CvsQ, name = "HSC_CvsQ", force = TRUE)
}

# chromVAR
proj <- addBgdPeaks(proj, force = TRUE)
proj <- addDeviationsMatrix(
  ArchRProj = proj,
  peakAnnotation = "HSC_CvsQ",
  force = TRUE
)
getAvailableMatrices(proj)

# save
saveArchRProject(ArchRProj = proj, load = FALSE)

Sys.Date()
sessionInfo()

