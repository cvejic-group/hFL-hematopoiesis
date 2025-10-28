#!/usr/bin/env Rscript

DOCNAME <- "archr_48FL"
dir.create(here::here("output", DOCNAME), showWarnings = FALSE)
# dirs
archr_dir <- here::here("output", DOCNAME)
gr_dir <- file.path(archr_dir, "PeakCalls")
bed_dir <- file.path(archr_dir, "PeakBED")

library(tidyverse)
library(ArchR)

gr_to_bed <- function(gr=NULL, out=NULL) {
  df <- data.frame(chr=seqnames(gr),
                   left=format(start(gr)-1, scientific=F, trim = T),
                   right=format(end(gr), scientific=F, trim = T),
                   names=".",
                   scores=gr$score,
                   strands=strand(gr),
                   distToGeneStart=gr$distToGeneStart,
                   nearestGene=gr$nearestGene,
                   peakType=gr$peakType,
                   distToTSS=gr$distToTSS,
                   nearestTSS=gr$nearestTSS,
                   GC=gr$GC) %>%
    mutate(names = paste0(chr, ":", left, "-", right))
  colnames(df)[1] <- "#chr"
  #return(df)
  write_tsv(df, file = out)
}



cells <- c("HSC", "GP", "Granulocyte",
           "MEMP-t", "MEMP", "MEP", "MEMP-Mast-Ery", "MEMP-Ery", "Early-Ery", "Late-Ery",
           "MEMP-MK", "MK", "MastP-t", "MastP", "Mast",
           "MDP", "Monocyte", "Kupffer", "cDC1", "cDC2", "pDC", "ASDC",
           "LMPP", "LP", "Cycling-LP", "PreProB", "ProB-1", "ProB-2", "Large-PreB", "Small-PreB", "IM-B",
           "NK", "ILCP", "T",
           "Hepatocyte", "Endothelia")
PCWs <- c(5:18)

if (FALSE) {
  # celltype
  dir.create(file.path(bed_dir, "celltype"), showWarnings = FALSE, recursive = TRUE)
  for (cell in cells) {
    cell2 <- str_replace_all(cell, "-", ".")
    cell2 <- case_when(
      cell == "MEMP-t" ~ "X49.MEMP",
      cell == "MDP" ~ "X35.MDP",
      cell == "MastP-t" ~ "X49.2.MastP",
      cell == "Cycling-LP" ~ "Cycling_LP",
      TRUE ~ cell2
    )
    print(cell2)
    peaks <- readRDS(file.path(gr_dir, "celltype", paste0(cell2, "-reproduciblePeaks.gr.rds")))
    out <- file.path(bed_dir, "celltype", paste0(cell, ".bed"))
    gr_to_bed(gr = peaks, out = out)
  }
}

# all peaks
ps <- readRDS("/work/home/project/20231127_DevM/gwas/sLDSC/atac_peaks/archr_peakSet.rds")
out <- file.path(bed_dir, "FL_all_peaks.bed")
gr_to_bed(gr = ps, out = out)


# PCW_celltype
dir.create(file.path(bed_dir, "PCW_celltype"), showWarnings = FALSE, recursive = TRUE)
for (cell in cells) {
  cell2 <- str_replace_all(cell, "-", ".")
  print(cell2)
  for (pcw in PCWs) {
    print(pcw)
    gr_file <- file.path(gr_dir, "PCW_celltype", paste0(paste0(cell2, "_PCW", pcw), "-reproduciblePeaks.gr.rds"))
    if (file.exists(gr_file)) {
      peaks <- readRDS(gr_file)
      out <- file.path(bed_dir, "PCW_celltype", paste0(paste0(cell, "_PCW", pcw), ".bed"))
      gr_to_bed(gr = peaks, out = out)
    }
  }
}

