#!/usr/bin/env Rscript

# run
# nohup Rscript analysis/archr.05.Peak2GeneLinks.R > analysis/archr.05.Peak2GeneLinks.log &

DOCNAME <- "archr_48FL"
dir.create(here::here("output", DOCNAME), showWarnings = FALSE)

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

proj <- addPeak2GeneLinks(
  ArchRProj = proj,
  reducedDims = "Harmony",
  useMatrix = "GeneExpressionMatrix",
  dimsToUse = 1:30,
  corCutOff = corCutOff
)

# save
saveArchRProject(ArchRProj = proj, load = FALSE)

R.utils::copyDirectory(here::here("output", "archr_48FL"), here::here("output", "archr_48FL.backup"))

Sys.Date()
sessionInfo()



