###################################
#                                 #
#          02.3 QC plots          #
#                                 #
###################################

## Project: Developmental multiome
## Date: 7.03.2022
## Author: Pernille
## Script name: 02.3.QC_plots.R

## Script description: 
## Cell filtering 

source("/work/home/Developmental_multiome/Analysis/Utils/0.Utils.R")

.libPaths("/work/home/Software/R/")

library(Seurat) 
library(Signac)
library(tidyverse)
library(tidyr)
library(dplyr)
library(gridExtra)

source("/work/DevM_analysis/01.annotation/02.cleandata/Utils.R")

##---------------------------------------------##
##----------------0. Load data-----------------##
##---------------------------------------------##

meta_cells <- read.table("/work/DevM_analysis/01.annotation/02.cleandata/data/FL_cell-info_filtering.txt", sep = "\t")
sam_info <- read.table(file = "/work/DevM_analysis/utils/sam_info.DevM.12.08.24.txt", sep = "\t", header = T)
sam_info <- sam_info[order(sam_info$PCW, decreasing =  F), ]

#For CB: 
#meta_cells <- read.table("/work/DevM_analysis/01.annotation/02.cleandata/data/CB_cell-info_filtering.txt", sep = "\t")
#sam_info <- read.table(file = "/work/DevM_analysis/utils/sam_info.DevM.20.08.24.txt", sep = "\t", header = T)
#sam_info <- sam_info[grep("CB", sam_info$libraryID), ]

meta_cells <- left_join(meta_cells, sam_info %>% select("libraryID", "sampleID", "donorID", "PCW", "Sex"),
                        by = c("libraryID", "sampleID", "donorID"))

##---------------------------------------------##
##-------------1. Barplot n cells--------------##
##---------------------------------------------##

## Barplot per donorID
df <- table(meta_cells$HighQualityCell,meta_cells$sampleID)
df <- df[,-grep("/", colnames(df))]
df_percent <- reshape2::melt(t(df)/colSums(df)*100)
colnames(df_percent) <- c("sample", "variable", "value")
df_percent$value <- round(df_percent$value, 2)
df <- reshape2::melt(t(df))
colnames(df) <- c("sample", "variable", "value")

p <- my_stacked_barplot(df, 
                   factor_levels_variables = c(0,1),
                   factor_level_samples = sam_info$sampleID, #in order of age
                   col = c("#bababa","#33a02c")) +
  aes_turn_x_text + xlab("donorID") + ylab("n cells") 
ggsave(plot = p, 
       filename = "/work/DevM_analysis/01.annotation/02.cleandata/plots/FL_Barplot_nHighQualityCells_perSampleID_orderPCW.pdf",
       height = 8, width = 24)

p <- my_stacked_barplot(df_percent, 
                        factor_levels_variables = c(0,1),
                        factor_level_samples = sam_info$sampleID, #in order of age
                        col = c("#bababa","#33a02c")) +
  aes_turn_x_text + xlab("donorID") + ylab("% cells") 

ggsave(plot = p, 
       filename = "/work/DevM_analysis/01.annotation/02.cleandata/plots/FL_Barplot_percentHighQualityCells_perSampleID_orderPCW.pdf",
       height = 8, width = 24)

## Barplot per PCW
df <- table(meta_cells$HighQualityCell,meta_cells$PCW)
df_percent <- reshape2::melt(t(df)/colSums(df)*100)
colnames(df_percent) <- c("sample", "variable", "value")
df_percent$value <- round(df_percent$value, 2)
df <- reshape2::melt(t(df))
colnames(df) <- c("sample", "variable", "value")

p <- my_stacked_barplot(df_percent, 
                        factor_levels_variables = c(0,1),
                        #factor_level_samples = sam_info$sampleID, #in order of age
                        col = c("#bababa","#33a02c")) +
  aes_turn_x_text + ylab("% cells") +
  scale_x_continuous(name ="PCW", 
                     breaks=6:18)
ggsave(plot = p, 
       filename = "/work/DevM_analysis/01.annotation/02.cleandata/plots/Barplot_percentHighQualityCells_perPCW.pdf",
       height = 7, width = 7)

p <- my_stacked_barplot(df, 
                        factor_levels_variables = c(0,1),
                        #factor_level_samples = sam_info$sampleID, #in order of age
                        col = c("#bababa","#33a02c")) +
  aes_turn_x_text + ylab("n cells") +
  scale_x_continuous(name ="PCW", 
                     breaks=6:18)
ggsave(plot = p, 
       filename = "/work/DevM_analysis/01.annotation/02.cleandata/plots/Barplot_nHighQualityCells_perPCW.pdf",
       height = 7, width = 7)

## Barplot per Sex
df <- table(meta_cells$HighQualityCell,meta_cells$Sex)
df <- df[, c("F", "M")]
df_percent <- reshape2::melt(t(df)/colSums(df)*100)
colnames(df_percent) <- c("sample", "variable", "value")
df_percent$value <- round(df_percent$value, 2)
df <- reshape2::melt(t(df))
colnames(df) <- c("sample", "variable", "value")

p <- my_stacked_barplot(df_percent, 
                        factor_levels_variables = c(0,1),
                        #factor_level_samples = sam_info$sampleID, #in order of age
                        col = c("#bababa","#33a02c")) +
  aes_turn_x_text + ylab("% cells") 
ggsave(plot = p, 
       filename = "/work/DevM_analysis/01.annotation/02.cleandata/plots/FL_Barplot_percentHighQualityCells_perSex.pdf",
       height = 3.6, width = 3)

p <- my_stacked_barplot(df, 
                        factor_levels_variables = c(0,1),
                        #factor_level_samples = sam_info$sampleID, #in order of age
                        col = c("#bababa","#33a02c")) +
  aes_turn_x_text + ylab("n cells")
ggsave(plot = p, 
       filename = "/work/DevM_analysis/01.annotation/02.cleandata/plots/FL_Barplot_nHighQualityCells_perSex.pdf",
       height = 3.6, width = 3)

## N cells per PCW and Sex
df <- meta_cells[-grep("/", meta_cells$donorID),]
df <- df %>% 
  filter(!Sex %in% c("Contaminated", "Unknown")) %>% 
  group_by(PCW, Sex, HighQualityCell) %>% 
  summarise(n = n())
df$HighQualityCell <- as.factor(df$HighQualityCell)

p <- ggplot(df) +
  geom_bar(aes(x = Sex, y = n, fill = HighQualityCell, col = Sex),
           position = "stack",
           stat = "identity") +
  facet_grid(~ PCW, switch = "x") +
  scale_color_manual(values = c("#e31a1c","#1f78b4")) +
  scale_fill_manual(values = c("#bababa","#33a02c")) + 
  ylab("n cells") + xlab(c("PCW\nSex"))+
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = NA, color = "white"),
        panel.spacing = unit(-.01,"cm")) +
  mashaGgplot2Theme
ggsave(plot = p, 
       filename = "/work/DevM_analysis/01.annotation/02.cleandata/plots/Barplot_nHighQualityCells_perPCW-Sex.pdf",
       height = 6.22, width = 12)

##---------------------------------------------##
##---------------2. Violin meta----------------##
##---------------------------------------------##
meta_cells$libraryID <- factor(meta_cells$libraryID, levels = unique(sam_info$libraryID)) #order by age
meta_cells$PCW <- as.factor(meta_cells$PCW)
P <- lapply(c("nCount_RNA", "nFeature_RNA", 
              "percent.mt", "percent.rb", 
              "nCount_peaks", "nFeature_peaks",
              "TSS.enrichment", "nucleosome_signal" ), 
            FUN = function(x){
              df <- meta_cells %>% rename("to_plot" = x)
              p <- ggplot(df %>% filter(HighQualityCell == 1), 
                          aes( x = libraryID, y = to_plot, fill = PCW)) + 
                scale_fill_manual(values = c("#620E3F",'#8e0152','#c51b7d','#de77ae','#f1b6da','#fde0ef','cornsilk',
                                             '#e6f5d0','#b8e186','#7fbc41','#4d9221','#276419', "#00441b", "#000A02")) + 
                geom_violin() + 
                mashaGgplot2Theme + theme(legend.position = "bottom") + aes_turn_x_text + ylab(x)
              return(p)})

ggsave(plot = marrangeGrob(P, ncol = 1, nrow = 1), 
       filename = "/work/DevM_analysis/01.annotation/02.cleandata/plots/FL_Vln_metaQC_HighQualityCells_perlibraryID_orderPCW.pdf",
       height = 8, width = 6)

##---------------------------------------------##
##-------------3. Umap per sample--------------##
##---------------------------------------------##

## Colored by meta data 
my_data <- list.files(path = "/work/DevM_analysis/01.annotation/02.cleandata/data/seu_obj", 
                      full.names = TRUE, recursive = F)

#For CB:
my_data <- my_data[grep("CB", my_data)]
P <- lapply(my_data, FUN = function(x){
  cat(x, "\n")
  seu <- readRDS(x)
  
  if(grepl("rna", x)){
    p <- FeaturePlot(seu, c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb"), order = T) 
      #scale_color_gradientn(colors = viridis::TurboPalette(50))
      #patchwork::plot_annotation(title = unique(seu$libraryID))
  }
  if(grepl("atac", x)){
    p <- FeaturePlot(seu, c("nCount_peaks", "nFeature_peaks", "TSS.enrichment", "nucleosome_signal"), order = T) 
      #scale_color_gradientn(colors = viridis::TurboPalette(50))
      #patchwork::plot_annotation(title = unique(seu$libraryID))
  }
  
  return(p)
})
names(P) <- str_remove(sapply(strsplit(my_data, "/"), '[',8), "_seu.rds")

pdf("/work/DevM_analysis/01.annotation/02.cleandata/plots/CB_UMAP_QCmetrics_perLibraryID.pdf", 
    width = 8.3, height = 6)
for(p in names(P)){
  print(P[[p]] + 
          patchwork::plot_annotation(title = paste(p, "-", sam_info %>% 
                                                     filter(libraryID == strsplit(p, "_")[[1]][1]) %>% 
                                                     select("PCW"), "PCW")))}
dev.off()


## Colored by donorID for pooled samples 
pooled_data <- my_data[grep("PCW", my_data)]

P <- lapply(pooled_data, FUN = function(x){
  cat(x, "\n")
  seu <- readRDS(x)
  p1 <- DimPlot(seu, group.by = "donorID", order = T)
  if(grepl("rna", x)){
    p2 <- DotPlot(seu, features = c("XIST",
                                    "TTTY15", "KDM5D", "TTTY14","ZFY", "PRKY", "GYG2P1", "PSMA6P1",
                                    "USP9Y", "TXLNGY", "UTY", "RPS4Y1", "DDX3Y", "EIF1AY"), group.by = "donorID")
    return(list(p1,p2))
  }else{
    return(p1)
  }
})
names(P) <- str_remove(sapply(strsplit(pooled_data, "/"), '[',8), "_seu.rds")

pdf("/work/DevM_analysis/01.annotation/02.cleandata/plots/UMAP_coloredDonorID_pooledSamples.pdf", 
    width = 8.3, height = 6)
for(p in names(P)){
  print(P[[p]][[1]] + 
          patchwork::plot_annotation(title = paste(p, "-", sam_info %>% 
                                                     filter(libraryID == strsplit(p, "_")[[1]][1]) %>% 
                                                     select("PCW"), "PCW")))}
dev.off()

pdf("/work/DevM_analysis/01.annotation/02.cleandata/plots/DotPlot_sexGenes_pooledSamples.pdf", 
    width = 8.3, height = 6)
for(p in names(P)){
  if(grepl("rna", p)){
    print(P[[p]][[2]] + 
            patchwork::plot_annotation(title = paste(p, "-", sam_info %>% 
                                                       filter(libraryID == strsplit(p, "_")[[1]][1]) %>% 
                                                       select("PCW"), "PCW")))} 
  }
dev.off()


## Colored by Markers 

## Colored by meta data 
my_data <- my_data[grep("rna", my_data)]

P <- lapply(my_data, FUN = function(x){
  cat(x, "\n")
  seu <- readRDS(x)
  
  p <- FeaturePlot(seu, features = c("THY1", "CD34", "MLLT3",
                                "GATA2", "MPO", "KIT",
                                "ITGA2B", "LYZ", "MME"), 
              order = T,  pt.size = 0.6)
  
  return(p)
})
names(P) <- str_remove(sapply(strsplit(my_data, "/"), '[',8), "_seu.rds")

pdf("/work/DevM_analysis/01.annotation/02.cleandata/plots/UMAP_fewQCmarkers_perLibraryID.pdf", 
    width = 8.3, height = 6)
for(p in names(P)){
  print(P[[p]] + 
          patchwork::plot_annotation(title = paste(p, "-", sam_info %>% 
                                                     filter(libraryID == strsplit(p, "_")[[1]][1]) %>% 
                                                     select("PCW"), "PCW")))}
dev.off()
