#!/usr/bin/env Rscript

# nohup Rscript code/02.filter_and_smooth.R > logs/02.filter_and_smooth.log &

DOCNAME <- "DNAm"
dir.create(here::here("output", DOCNAME), showWarnings = FALSE)

library(tidyverse)
library(data.table)
library(bsseq)
library(DelayedMatrixStats)
library(BiocParallel)
source(here::here("utils/DNAm.R"))

# load
bs_5mC <- qs2::qs_read(here::here("output", DOCNAME, "bs_5mC_raw.qs2"))
bs_5mC

# remove all CpGs that have no coverage in at least one sample
loci.idx <- which(DelayedMatrixStats::rowSums2(getCoverage(bs_5mC, type="Cov")==0) == 0)
bs_5mC_filt <- bs_5mC[loci.idx, ]
bs_5mC_filt

# smooth
bs_5mC_filt <- BSmooth(BSseq = bs_5mC_filt, BPPARAM = MulticoreParam(workers = 12), verbose = TRUE)

# sort by genomic coordinate
bs_5mC_filt <- sort(bs_5mC_filt)

# save
qs2::qs_save(bs_5mC_filt, file = here::here("output", DOCNAME, "bs_5mC_filt.qs2"))

# only autosomes
autosomes <- paste0("chr", 1:22)  # human autosomes
bs_5mC_auto <- keepSeqlevels(bs_5mC_filt, autosomes, pruning.mode = "coarse")
bs_5mC_auto

# sort by genomic coordinate
bs_5mC_auto <- sort(bs_5mC_auto)

# save
qs2::qs_save(bs_5mC_auto, file = here::here("output", DOCNAME, "bs_5mC_filt.autosomes.qs2"))

# only chr1 for test
bs_5mC_chr1 <- keepSeqlevels(bs_5mC_auto, "chr1", pruning.mode = "coarse")
bs_5mC_chr1
qs2::qs_save(bs_5mC_chr1, file = here::here("output", DOCNAME, "bs_5mC_filt.chr1.qs2"))


sessionInfo()
Sys.Date()

