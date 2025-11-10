#!/usr/bin/env Rscript

# nohup Rscript code/DNAm.genomic_elements.R > logs/DNAm.genomic_elements.log &

DOCNAME <- "DNAm"
dir.create(here::here("output", DOCNAME), showWarnings = FALSE)

library(tidyverse)
library(bsseq)
source(here::here("utils/DNAm.R"))

# bsseq
bs <- qs2::qs_read(here::here("output", DOCNAME, "bs_5mC_filt.autosomes.qs2"))

# promoter
gr_pro <- readRDS(here::here("data/cellranger-2020-A-2.0.0.promoter.rds"))
methy_pro <- getMeth(bs, regions = gr_pro, type = "raw", what = "perRegion")
dim(methy_pro)
saveRDS(methy_pro, file = here::here("output", DOCNAME, "5mC_promoter_methy.rds"))

# all ATAC peaks
gr_atac <- readRDS("/work/home/project/20231127_DevM/gwas/sLDSC/atac_peaks/archr_peakSet.rds")
methy_atac <- getMeth(bs, regions = gr_atac, type = "raw", what = "perRegion")
dim(methy_atac)
saveRDS(methy_atac, file = here::here("output", DOCNAME, "5mC_ATACpeak_methy.rds"))

# encode cCREs
gr_re <- ChIPseeker::readPeakFile("/work/DevM_analysis/data/refATAC/ENCODE_SCREEN/GRCh38-cCREs.bed")
methy_re <- getMeth(bs, regions = gr_re, type = "raw", what = "perRegion")
dim(methy_re)
saveRDS(methy_re, file = here::here("output", DOCNAME, "5mC_ENCODE_cCRE_methy.rds"))

# epimap k562 chromatin states
gr_cs <- ChIPseeker::readPeakFile("/work/DevM_analysis/data/refATAC/epimap/BSS00762_18_CALLS_segments.bed")
methy_cs <- getMeth(bs, regions = gr_cs, type = "raw", what = "perRegion")
dim(methy_cs)
saveRDS(methy_cs, file = here::here("output", DOCNAME, "5mC_EpiMap_K562_chromatinState_methy.rds"))

sessionInfo()
Sys.Date()






