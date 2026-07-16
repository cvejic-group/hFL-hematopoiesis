# ============================================================
# Title: HSC.dev_de_lmm_qScoreAdj.r
# Purpose:
#   Fit the linear mixed model to identify genes associated with
#   HSC quiescence while adjusting for quiescence and technical
#   covariates

# ============================================================
# Load packages
# ============================================================
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(tidyverse)
library(qs)
source("/work/Jupyterlab/Project/CB_new/06.FL_vs_CB_HSCs/LMM.R")

workdir <- '/work/Jupyterlab/Project/CellCycle/03.Disentanglement/'
datadir <- paste0(workdir,'data/')
plotdir <- paste0(workdir,'plot/')

# ============================================================
# Load data
# ============================================================
fl_hsc <- readRDS('/work/Jupyterlab/Project/CellCycle/02.Dynamic_CC/01.Identify_G0_HSCs/data/11.FLHSC_withG0label_byBMgset.rds')
qg <- readRDS('/work/Jupyterlab/Project/CellCycle/02.Dynamic_CC/01.Identify_G0_HSCs/data/quiescent_hsc_signatures.rds')
inputG <- readRDS('/work/DevM_analysis/data/HSC_dev_diff_genes/HSC_lmm_de_time_df.rds')$gene

# ============================================================
# Compute quiescence module score
# ============================================================
fl_hsc <- NormalizeData(fl_hsc)
fl_hsc <- ScaleData(fl_hsc,features = rownames(fl_hsc))
fl_hsc <- AddModuleScore(fl_hsc,features = list(qg$BM_LTvsST_up),
                         seed = 7,ctrl = 100, name = 'BM_LTvsST_up',search = TRUE)
fl_hsc

# ============================================================
# Preprocessing for LMM
# ============================================================

# filter genes
fl_hsc <- fl_hsc[,fl_hsc$PCW_new != 'Mixed']
fl_hsc <- fl_hsc[inputG,]

# tidy covariates
metadata <- fl_hsc@meta.data |> 
  dplyr::select(nCount_RNA,nFeature_RNA,percent.mt,percent.rb,
                PCW_new,libraryID,donorID,Sex,BM_LTvsST_up1,Batch)

metadata <- metadata |>
  rownames_to_column("barcode") |>
  mutate(
    logUMI = log(nCount_RNA),
    logFeature = log(nFeature_RNA),
    libraryID = factor(libraryID),
    donorID = factor(donorID),
    Sex = factor(Sex),
    PCW = as.numeric(scale(as.integer(as.character(PCW_new))))
  ) |>
  droplevels() |>
  column_to_rownames("barcode") |>
  as.data.frame()

# ============================================================
# Modeling
# ============================================================
table(metadata$PCW, useNA = "ifany")
str(metadata)

Y <- GetAssayData(fl_hsc, layer = "scale.data")
dim(Y)
stopifnot(all(colnames(Y) == rownames(metadata)))

res <- LFLMM(Y, metadata[,c("logUMI",'logFeature', "percent.mt","Batch", "libraryID", 
                            "donorID","Sex", "BM_LTvsST_up1","PCW")], ITRMAX=300)

saveRDS(res, file = paste0(datadir,"111.HSC_adjQ_lmm_res.rds"))

# de
de <- getBF(Y, res, "PCW", DE1 = NA)
names(de)

# de df
df_de <- data.frame(gene = rownames(de$beta), beta = de$beta[,1], ltsr = de$ltsr[,1]) |>
  arrange(-beta)

# ============================================================
# Save
# ============================================================
saveRDS(de, file = paste0(datadir,"111.HSC_adjQ_lmm_de_time.rds"))
saveRDS(df_de, file = paste0(datadir,"111.HSC_adjQ_lmm_de_time_df.rds"))
df_de_sig <- df_de[df_de$ltsr > 0.9,]
df_de_sig |>
  write_csv(file = paste0(datadir,"111.HSC_adjQ_lmm_de_time_ltsr0.9.csv"))

#
sessionInfo()
Sys.Date()




