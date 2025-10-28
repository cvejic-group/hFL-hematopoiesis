#!/usr/bin/env Rscript

# run
# nohup Rscript code/archr.09.subset_blood.R > logs/archr.09.subset_blood.log &

DOCNAME <- "archr_48FL.blood"
dir.create(here::here("output", DOCNAME), recursive = T, showWarnings = F)

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
proj <- loadArchRProject(path = here::here('output', "archr_48FL.slim"))
proj

# sub blood
idxCell <- BiocGenerics::which(proj$celltype %in% c("HSC", "MEMP-t", "MastP-t", "MDP", "LMPP", "LP", "Cycling-LP",
                                                    "GP", "Granulocyte",
                                                    "MEMP", "MEP", "MEMP-Mast-Ery", "MEMP-Ery", "Early-Ery", "Late-Ery",
                                                    "MEMP-MK", "MK", "MastP", "Mast",
                                                    "Monocyte", "Kupffer", "cDC1", "cDC2", "pDC", "ASDC",
                                                    "PreProB", "ProB-1", "ProB-2", "Large-PreB", "Small-PreB", "IM-B",
                                                    "NK", "ILCP", "T"))
cellsTarget <- proj$cellNames[idxCell]
proj_sub <- subsetArchRProject(ArchRProj = proj, cells = cellsTarget, dropCells = TRUE,
                               outputDirectory = here::here("output", DOCNAME), force = TRUE)
saveArchRProject(ArchRProj = proj_sub, load = FALSE)

