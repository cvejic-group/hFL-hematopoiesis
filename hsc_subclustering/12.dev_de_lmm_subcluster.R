#/usr/bin/env Rscript

# R432

library(tidyverse)
library(Seurat)
source("/work/home/RefData/SKM_ageing_atlas/DE_analysis/LMM.R")

DOCNAME <- "HSC.dev_de_lmm_subcluster"
dir.create(here::here("output", DOCNAME), showWarnings = FALSE)

# load mdata
mdata <- qs::qread(here::here("output", "HSC.metacell", "mdata.qs")) |>
  mutate(logUMI = log(nCount_RNA),
         logFeature = log(nFeature_RNA),
         PCWsca = as.numeric(scale(as.integer(as.character(PCW))))) |>
  droplevels() |>
  as.data.frame()
str(mdata)

# gene mat
Y <- qs::qread(here::here("output", "HSC.metacell", "gm_scale.qs"))
dim(Y)
stopifnot(all(colnames(Y) == rownames(mdata)))

# de for sub
for (i in levels(mdata$celltype)) {
  print(i)
  # sub
  mdata_sub <- mdata |>
    filter(celltype == i) |>
    droplevels()
  Y_sub <- Y[,rownames(mdata_sub)]

  # linear mixed model
  res <- LFLMM(Y_sub, mdata_sub[,c("logUMI", "logFeature", "percent.mt", "libraryID", "donorID",
                           "Sex", "Batch", "PCWsca")], ITRMAX=300)
  saveRDS(res, file = here::here("output", DOCNAME, "HSC_lmm_res.rds"))
  # de
  de <- getBF(Y_sub, res, "PCWsca", DE1 = NA)
  names(de)
  # de df
  df_de <- data.frame(gene = rownames(de$beta), beta = de$beta[,1], ltsr = de$ltsr[,1]) |>
    arrange(-beta)
  saveRDS(de, file = here::here("output", DOCNAME, paste0(i, "_lmm_de_time.rds")))
  saveRDS(df_de, file = here::here("output", DOCNAME, paste0(i, "_lmm_de_time_df.rds")))
  df_de_sig <- df_de[df_de$ltsr > 0.9,]
  df_de_sig |>
    write_csv(file = here::here("output", DOCNAME, paste0(i, "_lmm_de_time_ltsr0.9.csv")))
}

sessionInfo()
Sys.Date()
