#### Set Up ####
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

##### load  data #####
peak <- readRDS(paste0(datadir,'21.Peak_clusterbyphase_kmeans.rds'))
peak_gr <- lapply(peak,function(p){
  peakbed <- GRanges(p)
  peakbed
})
names(peak_gr) <- names(peak)

for(name in names(peak_gr)){
  peak_gr[[name]]$phase <- name
}
peak_grs <- unlist(GRangesList(peak_gr), use.names = FALSE)

#### monalisa motif enrichment ####
# Bin genomic regions to correct compositional bias
bins <- factor(peak_grs$phase,levels = names(peak_gr))
table(bins)


levels(bins)

# use latest motif annotation
jaspar <- JASPAR2024()
sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), db(jaspar))
pwms <- TFBSTools::getMatrixSet(sq24, list(matrixtype = "PWM", tax_group = "vertebrates", collection = "CORE",
                                           species='Homo sapiens'))

# check gc
summary(width(peak_grs))
seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, peak_grs)
plotBinDiagnostics(seqs = seqs, bins = bins, aspect = "GCfrac")
plotBinDiagnostics(seqs = seqs, bins = bins, aspect = "dinucfreq")

# binned motif
se1 <- calcBinnedMotifEnrR(seqs = seqs, bins = bins, pwmL = pwms,
                           verbose = TRUE)
#se1 <- readRDS(paste0(datadir,'24.Peak_clusteredbyphase_monalisa.rds'))
sel1 <- apply(assay(se1, "negLog10Padj"), 1, 
              function(x) max(abs(x), 0, na.rm = TRUE)) > 4
sum(sel1)
seSel1 <- se1[sel1, ]
Cairo::CairoPDF(paste0(plotdir,'24.Peak_clusteredbyphase_monalisa.pdf'), 
                width = 10, height = 13,family = 'Arial')
plotMotifHeatmaps(x = seSel1, which.plots = c("log2enr", "negLog10Padj"), 
                  width = 2.0, cluster = TRUE, maxEnr = 2, maxSig = 10, show_seqlogo = TRUE,
                  show_motif_GC = TRUE,show_bin_legend = TRUE)
dev.off()

# save
saveRDS(se1,paste0(datadir,'24.Peak_clusteredbyphase_monalisa.rds'))


