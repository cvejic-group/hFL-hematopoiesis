#!/usr/bin/env Rscript

# conda activate R412
# nohup Rscript 01.prep_diff_peak.R > 01.prep_diff_peak.log &

# avoid scientific notation
options(scipen = 999)

# set up
work_dir <- "~/project/20231127_DevM/gwas/sLDSC"
out_dir <- file.path(work_dir, "atac_peaks")


library(tidyverse)
library(ArchR)

set.seed(1)
# default number of threads
addArchRThreads(threads = 16)
# add a reference genome annotation for ArchR
addArchRGenome("hg38")
# chr prefix to exclude KI/GL scaffolds
addArchRChrPrefix(chrPrefix = TRUE)
# cor cutoff
corCutOff = 0.5

# use this backup for the peakSet
# as the peakSet of new archr proj is from "per celltype per PCW"
proj_dir <- "/work/DevM_analysis/data/ArchR_backup/archr_48FL"
proj <- loadArchRProject(path = proj_dir)
proj

# cluster-specific peaks
markersPeaks <- readRDS(file.path(proj_dir, "MarkerFeatures", "markersPeaks.rds"))
diffPeaks_FDR5_LFC1 <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 1")

save_diff_peak <- function(diffPeaks = NULL, save_dir = NULL, padding = c(-1000, 1000)) {
  # use new names
  col_order <- c("HSC", "GP", "Granulocyte",
                 "49-MEMP", "MEMP", "MEP", "MEMP-Mast-Ery", "MEMP-Ery", "Early-Ery", "Late-Ery",
                 "MEMP-MK", "MK", "49:2-MastP", "MastP", "Mast", "35-MDP",
                 "Monocyte", "Kupffer", "cDC1", "cDC2", "pDC", "ASDC",
                 "LMPP", "LP", "Cycling_LP", "PreProB", "ProB-1", "ProB-2",
                 "Large-PreB", "Small-PreB", "IM-B",
                 "NK", "ILCP", "T",
                 "Hepatocyte", "Endothelia")
  col_name_new <- c("HSC", "GP", "Granulocyte",
                    "MEMP-t", "MEMP", "MEP", "MEMP-Mast-Ery", "MEMP-Ery", "Early-Ery", "Late-Ery",
                    "MEMP-MK", "MK", "MastP-t", "MastP", "Mast", "MDP",
                    "Monocyte", "Kupffer", "cDC1", "cDC2", "pDC", "ASDC",
                    "LMPP", "LP", "Cycling-LP", "PreProB", "ProB-1", "ProB-2",
                    "Large-PreB", "Small-PreB", "IM-B",
                    "NK", "ILCP", "T",
                    "Hepatocyte", "Endothelia")
  for (i in 1:length(col_order)){
    x = col_order[i]
    name = col_name_new[i]
    #print(name)
    differtest <- diffPeaks[[x]]
    differtest <- as.data.frame(differtest) %>%
      filter(seqnames %in% paste0("chr", 1:22)) %>%
      dplyr::select(seqnames, start, end) %>%
      mutate(start = as.integer(start - padding[1]),
             end =  as.integer(end + padding[2])) %>%
      mutate(start = pmax(start, 0))
    out_path = file.path(out_dir, save_dir, name)
    dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
    out_name = file.path(out_path, paste0(name, "_peak.hg38.bed"))
    write_tsv(differtest, file=out_name, col_names = FALSE)
  }
}

# FDR5_LFC1 padding 1kb
save_diff_peak(diffPeaks = diffPeaks_FDR5_LFC1, save_dir = "FDR5_LFC1_padding_1k", padding = c(-1000, 1000))

# all peaks as background
gr <- proj@peakSet
length(proj@peakSet)

# ctrl peaks to bed
ctrl_to_bed <- function(gr=NULL, save_dir=NULL, padding = c(-1000, 1000)) {
  df <- data.frame(chr=seqnames(gr),
                   left=start(gr),
                   right=end(gr)) %>%
    filter(chr %in% paste0("chr", 1:22)) %>%
    mutate(left = left - padding[1],
           right = right + padding[2]) %>%
    mutate(left = pmax(left, 0))
  out_path <-file.path(out_dir, save_dir, "CTRL")
  dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
  df %>%
    write_tsv(file=file.path(out_path, "CTRL_peak.hg38.bed"), col_names = FALSE)
}

# save all peaks
saveRDS(gr, file = file.path(out_dir, "archr_peakSet.rds"))

# to different dir
ctrl_to_bed(gr = gr, save_dir = "FDR5_LFC1_padding_1k", padding = c(-1000, 1000))

# session info
devtools::session_info()
Sys.Date()
