#!/usr/bin/env Rscript

# run
# nohup Rscript code/archr.getGroupSE.R > logs/archr.getGroupSE.log &

DOCNAME <- "archr_48FL.blood"
# dirs
archr_dir <- here::here("output", DOCNAME)
se_dir <- file.path(archr_dir, "GroupSE")
dir.create(se_dir, showWarnings = FALSE, recursive = TRUE)

library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyverse)
library(ArchR)
#reticulate::use_python("/work/home/software/anaconda3/bin/python")

# set up
set.seed(1)
# default number of threads
addArchRThreads(threads = 4)
# add a reference genome annotation for ArchR
addArchRGenome("hg38")
# chr prefix to exclude KI/GL scaffolds
addArchRChrPrefix(chrPrefix = TRUE)
# cor cutoff
corCutOff = 0.5

# load
proj <- loadArchRProject(path = here::here('output', DOCNAME))
proj
#getAvailableMatrices(proj)

# group by sample
se <- getGroupSE(ArchRProj = proj, useMatrix = "PeakMatrix", groupBy = "sampleID", divideN = FALSE, scaleTo = NULL)
saveRDS(se, file = file.path(se_dir, "PeakMatrix.PB_sampleID.rds"))

# get fraction of open of each peak
arrow_files <- getArrowFiles(proj)
peak_counts_list <- list()
# Loop through Arrow files
for (af in arrow_files) {
  # Read the sparse PeakMatrix for this Arrow file
  peak_data <- getMatrixFromArrow(
    ArrowFile = af,
    useMatrix = "PeakMatrix",
    binarize = TRUE      # Only care about presence/absence
  )
  peak_counts <- Matrix::rowSums(assay(peak_data) > 0)
  peak_counts_list[[af]] <- peak_counts
}

# sum the counts for each peak (names will align)
combined_peak_counts <- Reduce("+", peak_counts_list)

# all cells
total_cells <- 304786
peak_fraction <- combined_peak_counts / total_cells
peak_name <- paste0(as.character(seqnames(proj@peakSet)), '-',
                    as.character(start(proj@peakSet)), '-', as.character(end(proj@peakSet)))
write_csv(
  data.frame(Peak=peak_name, Fraction=peak_fraction),
  file = file.path(se_dir, "PeakFractionInCells.tsv")
)

sessionInfo()
Sys.Date()
