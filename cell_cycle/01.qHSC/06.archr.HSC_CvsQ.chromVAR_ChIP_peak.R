#!/usr/bin/env Rscript

# run
# nohup Rscript code/archr.HSC_CvsQ.chromVAR_ChIP_peak.R > logs/archr.HSC_CvsQ.chromVAR_ChIP_peak.log &

DOCNAME <- "archr_48FL"

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

data_dir <- "archr_48FL.eachCell"
cell <- "HSC"

# load
proj_dir <- here::here('output', data_dir, cell)
proj <- loadArchRProject(path = proj_dir)
proj
getAvailableMatrices(proj)

# ChIP peaks
ChIP_peak <- c(
  MEP_CTCF = "/work/home/project/20231127_DevM/HSC.CvsQ.peaks/MEP_ChIP/CTCF.bed",
  MEP_ERG = "/work/home/project/20231127_DevM/HSC.CvsQ.peaks/MEP_ChIP/ERG.bed",
  MEP_FLI1 = "/work/home/project/20231127_DevM/HSC.CvsQ.peaks/MEP_ChIP/FLI1.bed",
  MEP_GATA2 = "/work/home/project/20231127_DevM/HSC.CvsQ.peaks/MEP_ChIP/GATA2.bed",
  MEP_LMO2 = "/work/home/project/20231127_DevM/HSC.CvsQ.peaks/MEP_ChIP/LMO2.bed",
  MEP_LYL1 = "/work/home/project/20231127_DevM/HSC.CvsQ.peaks/MEP_ChIP/LYL1.bed",
  MEP_Pol2 = "/work/home/project/20231127_DevM/HSC.CvsQ.peaks/MEP_ChIP/Pol2.bed",
  MEP_PU1 = "/work/home/project/20231127_DevM/HSC.CvsQ.peaks/MEP_ChIP/PU1.bed",
  MEP_RUNX1 = "/work/home/project/20231127_DevM/HSC.CvsQ.peaks/MEP_ChIP/RUNX1.bed",
  MEP_STAG2 = "/work/home/project/20231127_DevM/HSC.CvsQ.peaks/MEP_ChIP/STAG2.bed",
  MEP_TAL1 = "/work/home/project/20231127_DevM/HSC.CvsQ.peaks/MEP_ChIP/TAL1.bed",
  GMP_CTCF = "/work/home/project/20231127_DevM/HSC.CvsQ.peaks/GMP_ChIP/CTCF.bed",
  GMP_ERG = "/work/home/project/20231127_DevM/HSC.CvsQ.peaks/GMP_ChIP/ERG.bed",
  GMP_FLI1 = "/work/home/project/20231127_DevM/HSC.CvsQ.peaks/GMP_ChIP/FLI1.bed",
  GMP_GATA2 = "/work/home/project/20231127_DevM/HSC.CvsQ.peaks/GMP_ChIP/GATA2.bed",
  GMP_LMO2 = "/work/home/project/20231127_DevM/HSC.CvsQ.peaks/GMP_ChIP/LMO2.bed",
  GMP_LYL1 = "/work/home/project/20231127_DevM/HSC.CvsQ.peaks/GMP_ChIP/LYL1.bed",
  GMP_Pol2 = "/work/home/project/20231127_DevM/HSC.CvsQ.peaks/GMP_ChIP/Pol2.bed",
  GMP_PU1 = "/work/home/project/20231127_DevM/HSC.CvsQ.peaks/GMP_ChIP/PU1.bed",
  GMP_RUNX1 = "/work/home/project/20231127_DevM/HSC.CvsQ.peaks/GMP_ChIP/RUNX1.bed",
  GMP_STAG2 = "/work/home/project/20231127_DevM/HSC.CvsQ.peaks/GMP_ChIP/STAG2.bed",
  GMP_TAL1 = "/work/home/project/20231127_DevM/HSC.CvsQ.peaks/GMP_ChIP/TAL1.bed"
)

if("HSPC_ChIP" %ni% names(proj@peakAnnotation)){
  proj <- addPeakAnnotations(ArchRProj = proj, regions = ChIP_peak, name = "HSPC_ChIP", force = TRUE)
}

# chromVAR
proj <- addBgdPeaks(proj, force = TRUE)
proj <- addDeviationsMatrix(
  ArchRProj = proj,
  peakAnnotation = "HSPC_ChIP",
  force = TRUE
)
getAvailableMatrices(proj)

# save
saveArchRProject(ArchRProj = proj, load = FALSE)

Sys.Date()
sessionInfo()

