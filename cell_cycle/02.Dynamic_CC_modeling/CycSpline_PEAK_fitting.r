#### Set Up ####
library(mgcv)
library(Matrix)
library(BiocParallel)
library(tidyverse)
library(progressr)
library(doFuture)  
library(doRNG)    
library(future)
library(ggplot2)

workdir <- '/work/Jupyterlab/Project/CellCycle/02.Dynamic_CC/06.GAM_CycSpline/'
datadir <- paste0(workdir,'data/')
plotdir <- paste0(workdir,'plots/')

#### Load Data ####
peak_meta <- readRDS('/work/Jupyterlab/Project/CellCycle/02.Dynamic_CC/04.Dynamic_CC_ATAC/data/01.PeakMatrix_metacell.rds')
cc_meta <- read.csv('/work/Jupyterlab/Project/CellCycle/02.Dynamic_CC/03.Dynamic_CC_GEX/data/01.GEXMatrix_metacell_metaid.csv')
cc_meta <- cc_meta |> dplyr::distinct(meta_id,.keep_all = TRUE) |> 
  dplyr::select(meta_id,meta_coor,PCW_new) |>
  column_to_rownames('meta_id')
cc_meta <- cc_meta |> arrange(meta_coor)
peak_meta <- peak_meta[,rownames(cc_meta)]
all(colnames(peak_meta) == rownames(cc_meta))

#### fit fourier ####
source(paste0(workdir,'00.GAM_CycSpline_model.r'))

libsize <- Matrix::colSums(peak_meta)
res <- fit_all(counts_matrix = peak_meta,
               pseudotime = cc_meta$meta_coor,
               pcw = cc_meta$PCW_new,
               libsize = libsize,
               knots = 10,
               k = 40,
               workers = 8,
               method = "ML"
)

saveRDS(res,paste0(datadir,'20.Peak_fit_CycSpline.rds'))






