#!/usr/bin/env Rscript

# run
# nohup Rscript code/archr.getGroupBW.R > logs/archr.getGroupBW.log &

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

# by Library
#getGroupBW(ArchRProj = proj, groupBy = "libraryID", normMethod = "ReadsInTSS", tileSize = 50)

# by Donor
#getGroupBW(ArchRProj = proj, groupBy = "donorID", normMethod = "ReadsInTSS", tileSize = 50)

# by cell type
#getGroupBW(ArchRProj = proj, groupBy = "celltype", normMethod = "ReadsInTSS", tileSize = 50)

# by PCW_celltype
getGroupBW(ArchRProj = proj, groupBy = "PCW_celltype", normMethod = "ReadsInTSS", tileSize = 50)


