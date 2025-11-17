#### Set Up ####
library(gridExtra)
library(grid)
library(clusterProfiler)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(ReactomePA)
library(annotables)
library(org.Hs.eg.db)
library(patchwork)
library(ggplot2)
library(Cairo)
library(ggsci)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
set.seed(777)

workdir <- '/work/Jupyterlab/Project/CellCycle/02.Dynamic_CC/06.GAM_CycSpline/'
datadir <- paste0(workdir,'data/')
plotdir <- paste0(workdir,'plots/')

#### Load data ####
peak <- readRDS(paste0(datadir,'21.Peak_clusterbyphase_kmeans.rds'))
peak_gr <- lapply(peak,function(p){
  peakbed <- GRanges(p)
  peakbed
})
names(peak_gr) <- names(peak)

#### annotation ####
peak_anno <- list()
for(phase in names(peak_gr)){
  print(phase)
  gr <- peak_gr[[phase]]
  peak_anno[phase] <- annotatePeak(gr, tssRegion=c(-2000, 2000),
                                   TxDb=txdb, annoDb="org.Hs.eg.db")
}

Cairo::CairoPDF(paste0(plotdir,'22.plotDistToTSS.pdf'), 
                width = 10, height = 5,family='Arial')
plotDistToTSS(peak_anno) + scale_fill_jco() + theme(panel.grid = element_blank())
dev.off()

# save
saveRDS(peak_anno,paste0(datadir,'22.Peak_ClusteredByPhase_annotation.rds'))


####  plot #### 
peak_anno <- readRDS(paste0(datadir,'22.Peak_ClusteredByPhase_annotation.rds'))

d4p_list <- list()
for(phase in names(peak_anno)){
  anno <- peak_anno[[phase]]
  df<- getAnnoStat(anno) 
  df <- df |>
    mutate(Count = floor(Frequency * length(peak_gr[[phase]]) / 100),
           status = rep(phase,dim(df)[1]))
  d4p_list[[phase]] <- df
}

d4p <- do.call(rbind, d4p_list)

cols <- pal_flatui(palette = c("default"), alpha = 1)(10)
cols_val <- c(
  'Promoter (<=1kb)'      = cols[1],
  'Promoter (1-2kb)'      = cols[2],
  "5' UTR"                = cols[3],
  "3' UTR"                = cols[4],
  '1st Exon'              = cols[5],
  'Other Exon'            = cols[6],
  '1st Intron'            = cols[7],
  'Other Intron'          = cols[8],
  'Downstream (<=300)'    = cols[9],
  'Distal Intergenic'     = cols[10]
)
lvl <- names(cols_val)
d4p$Feature <- factor(d4p$Feature, levels = rev(levels(d4p$Feature)))
d4p$status <- factor(d4p$status, levels=names(peak_gr))

Cairo::CairoPDF(paste0(plotdir,'22.CCphase_peaknumber_barplot.pdf'), 
                width = 15, height = 7,family='Arial')
ggplot(d4p,aes(x=status,y=Count,fill=Feature)) + 
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  scale_fill_manual(values=cols_val) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA,color = 'black'),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)) +
  scale_x_discrete(drop = FALSE) +
  scale_x_discrete(breaks=names(peak_gr),
                   labels=c('mg1(180)','g1(1587)','g1s(2357)',
                            's(1629)','sg2(2682)','g2m(2381)')) +
  scale_y_continuous(breaks = c(500,1000,1500,2000,2500)) +
  xlab('') + ylab('Count')
dev.off()

# proportion 
Cairo::CairoPDF(paste0(plotdir,'22.CCphase_peakproportion_barplot.pdf'), 
                width = 15, height = 7,family='Arial')
ggplot(d4p,aes(x=status,y=Frequency,fill=Feature)) + 
  geom_bar(stat = "identity", position = "stack") +
  coord_flip() +
  scale_fill_manual(values=cols_val) + 
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA,color = 'black'),
        axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)) +
  scale_x_discrete(drop = FALSE)+
  scale_x_discrete(breaks=names(peak_gr),
                   labels=c('mg1(180)','g1(1587)','g1s(2357)',
                            's(1629)','sg2(2682)','g2m(2381)')) +
  xlab('') + ylab('Proportion')
dev.off()

