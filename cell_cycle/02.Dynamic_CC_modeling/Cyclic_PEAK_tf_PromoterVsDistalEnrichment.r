# ============================================================
# Title: Cyclic_PEAK_tf_PromoterVsDistalEnrichment.r
# Purpose:
#   Motif enrichment analysis on cyclic peaks separated by 
#   promoter and distal regions in cycling HSCs across cell
#   cycle phases

# ============================================================
# Load packages
# ============================================================
library(dplyr)
library(stringr)
library(ggplot2)
library(GenomicRanges)
library(SummarizedExperiment)
library(JASPAR2024)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(monaLisa)
library(ComplexHeatmap)
library(tidyverse)

workdir <- '/work/Jupyterlab/Project/CellCycle/02.Dynamic_CC/06.GAM_CycSpline/'
datadir <- paste0(workdir,'data/')
plotdir <- paste0(workdir,'plots/')

# ============================================================
# Jaspar 2024
# ============================================================
jaspar <- JASPAR2024()
sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), db(jaspar))
pwms <- TFBSTools::getMatrixSet(sq24, list(matrixtype = "PWM", tax_group = "vertebrates", collection = "CORE",
                                           species='Homo sapiens'))

# ============================================================
# monalisa
# ============================================================
# Get proximal and distal peak sets
peak_anno <- readRDS(paste0(datadir,'22.Peak_ClusteredByPhase_annotation.rds'))
proximal <- list()
distal <- list()
for(phase in names(peak_anno)){
  df <- as.data.frame(peak_anno[[phase]]@anno)
  p_peaks <- df |> filter(annotation %in% c('Promoter (<=1kb)'
                                            #'Promoter (1-2kb)'
  )) |>
    dplyr::select(seqnames,start,end) |>
    mutate(phase=phase)
  proximal[[phase]] <- p_peaks
  
  d_peaks <- df |> filter(annotation %in% c(#'1st Intron',
    #'Other Intron',
    'Distal Intergenic')) |>
    dplyr::select(seqnames,start,end) |>
    mutate(phase=phase)
  distal[[phase]] <- d_peaks
}
lapply(proximal, function(p){dim(p)})
lapply(distal, function(p){dim(p)})


proximal_gr <- unlist(GRangesList(proximal), use.names = FALSE)
distal_gr <- unlist(GRangesList(distal), use.names = FALSE)

bins_promoter <- factor(proximal_gr$phase,levels = names(proximal))
bins_distal <- factor(distal_gr$phase,levels = names(distal))
table(bins_distal)
table(bins_promoter)

seqs_promoter <- getSeq(BSgenome.Hsapiens.UCSC.hg38, proximal_gr)
seqs_distal <- getSeq(BSgenome.Hsapiens.UCSC.hg38, distal_gr)

se_promoter <- calcBinnedMotifEnrR(seqs = seqs_promoter, bins = bins_promoter, pwmL = pwms,
                                   verbose = TRUE)
se_distal <- calcBinnedMotifEnrR(seqs = seqs_distal, bins = bins_distal, pwmL = pwms,
                                 verbose = TRUE)
saveRDS(se_promoter, paste0(datadir,'24.Peak_promoterwithin1kb_monalisa.rds'))
saveRDS(se_distal, paste0(datadir,'24.Peak_distalintergenic_monalisa.rds'))

