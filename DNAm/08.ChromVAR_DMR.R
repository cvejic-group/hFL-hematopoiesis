#!/usr/bin/env Rscript

# run
# nohup Rscript code/08.ChromVAR_DMR.R

DOCNAME <- "archr_48FL.8to18"

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

# DMRs
DSS_DMR <- c(
  hyperDMR = "/work/home/project/20231127_DevM/DNAm/5mC_DSS_DMR_p005_Up.bed",
  hypoDMR = "/work/home/project/20231127_DevM/DNAm/5mC_DSS_DMR_p005_Dn.bed"
)

if("DSSDMR" %ni% names(proj@peakAnnotation)){
  proj <- addPeakAnnotations(ArchRProj = proj, regions = DSS_DMR, name = "DSSDMR", force = TRUE)
}

# chromVAR
proj <- addBgdPeaks(proj, force = TRUE)
proj <- addDeviationsMatrix(
  ArchRProj = proj,
  peakAnnotation = "DSSDMR",
  force = TRUE
)
getAvailableMatrices(proj)

# save
saveArchRProject(ArchRProj = proj, load = FALSE)

Sys.Date()
sessionInfo()

