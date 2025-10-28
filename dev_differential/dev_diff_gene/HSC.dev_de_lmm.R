#/usr/bin/env Rscript

# R432
# cd ~/project/20231127_DevM/devm_r432
# nohup Rscript code/HSC.dev_de_lmm.R > logs/HSC.dev_de_lmm.log &

library(tidyverse)
library(Seurat)
source("/work/home/software/SKM_ageing_atlas/DE_analysis/LMM.R")

DOCNAME <- "HSC.dev_de_lmm"
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

# linear mixed model
res <- LFLMM(Y, mdata[,c("logUMI", "logFeature", "percent.mt", "libraryID", "donorID",
                         "Sex", "Batch", "PCWsca")], ITRMAX=300)
saveRDS(res, file = here::here("output", DOCNAME, "HSC_lmm_res.rds"))
# de
de <- getBF(Y, res, "PCWsca", DE1 = NA)
names(de)
# de df
df_de <- data.frame(gene = rownames(de$beta), beta = de$beta[,1], ltsr = de$ltsr[,1]) |>
  arrange(-beta)
saveRDS(de, file = here::here("output", DOCNAME, "HSC_lmm_de_time.rds"))
saveRDS(df_de, file = here::here("output", DOCNAME, "HSC_lmm_de_time_df.rds"))
df_de_sig <- df_de[df_de$ltsr > 0.9,]
df_de_sig |>
  write_csv(file = here::here("output", DOCNAME, "HSC_lmm_de_time_ltsr0.9.csv"))


sessionInfo()
Sys.Date()
