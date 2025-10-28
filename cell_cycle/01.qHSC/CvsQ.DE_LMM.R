#!/usr/bin/env Rscript

# nohup Rscript code/HSC.CvsQ.DE_lmm.R > logs/HSC.CvsQ.DE_lmm.log &

DOCNAME <- "HSC.CvsQ.DE"
dir.create(here::here("output", DOCNAME), showWarnings = FALSE)

library(tidyverse)
library(Seurat)
source("/work/home/software/SKM_ageing_atlas/DE_analysis/LMM.R")

srat <- readRDS("/work/DevM_analysis/05.cellcycle/00.Identify_G0_HSCs/data/FLHSC_withG0label_byBMgset.rds")
table(srat$G0_label, useNA = "ifany")

# mdata
mdata <- srat@meta.data |>
  rownames_to_column("barcode") |>
  dplyr::select(-celltype, -starts_with("RNA_snn")) |>
  mutate(Q = ifelse(G0_label == "G0", "G0", "Other")) |>
  mutate(
    logUMI = log(nCount_RNA),
    logFeature = log(nFeature_RNA),
    Q = factor(Q, levels = c("Other", "G0")),
    donorID = factor(donorID),
    Sex = factor(Sex),
    Batch = factor(Batch)
  ) |>
  droplevels() |>
  column_to_rownames("barcode") |>
  as.data.frame()
dim(mdata)
table(mdata$Q, useNA = "ifany")
str(mdata)

# exp data (1 count 1% single-cell)
# scale Y
flag_rna <- qs::qread(here::here("output", "HSC.metacell", "hsc_gene.flag_1_1.qs"))
srat <- ScaleData(srat, features = rownames(srat)[flag_rna])
Y <- GetAssayData(srat, layer = "scale.data")
dim(Y)
stopifnot(all(colnames(Y) == rownames(mdata)))

# linear mixed model
res <- LFLMM(Y, mdata[,c("logUMI", "logFeature", "percent.mt", "libraryID", "donorID",
                         "Sex", "Batch", "Q")], ITRMAX=300)
saveRDS(res, file = here::here("output", DOCNAME, "G0_de.by_lmm.res.rds"))

# de
de <- getBF(Y, res, "Q", DE1 = 1)
names(de)
head(de$beta)
head(de$ltsr)

# de df
df_de <- data.frame(gene = rownames(de$beta), beta = de$beta[,2], ltsr = de$ltsr[,1]) |>
  arrange(-beta)
saveRDS(de, file = here::here("output", DOCNAME, "G0_de.by_lmm.res_de.rds"))
saveRDS(df_de, file = here::here("output", DOCNAME, "G0_de.by_lmm.res_df.rds"))
df_de_sig <- df_de[df_de$ltsr > 0.9,]
df_de_sig |>
  write_csv(file = here::here("output", DOCNAME, "G0_de.by_lmm.res_ltsr0.9.csv"))

# num
table(ifelse(df_de_sig$beta > 0, "Up", "Dn"))

sessionInfo()
Sys.Date()

