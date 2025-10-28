#/usr/bin/env Rscript

# R432
# cd /work/home/project/20231127_DevM/devm_r432
# nohup Rscript code/HSC.dev_da_lmm_mc_scale.R > logs/HSC.dev_da_lmm_mc_scale.log &

library(tidyverse)
library(Seurat)
library(Matrix)
source("/work/home/software/SKM_ageing_atlas/DE_analysis/LMM.R")

DOCNAME <- "HSC.dev_da_lmm_mc_scale"
dir.create(here::here("output", DOCNAME), showWarnings = FALSE)

# load mdata
mdata <- qs::qread(here::here("output", "HSC.metacell", "mc_mdata.qs")) |>
  mutate(PCWsca = as.numeric(scale(as.integer(as.character(PCW))))) |>
  droplevels() |>
  as.data.frame()
str(mdata)

# peak mat
Y_scale <- qs::qread(here::here("output", "HSC.metacell", "mc_pm_scale.qs"))
stopifnot(all(colnames(Y_scale) == mdata$mc_group))

# linear mixed model
res <- LFLMM(Y_scale, mdata[,c("m_frag", "libraryID", "donorID",
                               "Sex", "Batch", "PCWsca")], ITRMAX=300)
saveRDS(res, file = here::here("output", DOCNAME, "HSC_lmm_res.rds"))

# de with scale.data
if (TRUE) {
  # de
  de <- getBF(Y_scale, res, "PCWsca", DE1 = NA)
  names(de)
  # de df
  df_de <- data.frame(gene = rownames(de$beta), beta = de$beta[,1], ltsr = de$ltsr[,1]) |>
    arrange(-beta)
  saveRDS(de, file = here::here("output", DOCNAME, "HSC_dev_da_lmm.rds"))
  saveRDS(df_de, file = here::here("output", DOCNAME, "HSC_dev_da_lmm_df.rds"))
  # export peaks
  df_de |>
    mutate(name = gene) |>
    separate("gene", c("chr", "left", "right"), sep = "-") |>
    dplyr::select(chr, left, right, name, beta, ltsr) |>
    write_tsv(file = here::here("output", DOCNAME, "HSC_dev_da_lmm.peak_tested.bed"), col_names = FALSE)
  df_de |>
    mutate(name = gene) |>
    filter(ltsr > 0.9, beta > 0) |>
    separate("gene", c("chr", "left", "right"), sep = "-") |>
    dplyr::select(chr, left, right, name, beta) |>
    write_tsv(file = here::here("output", DOCNAME, "HSC_dev_da_lmm.peak_up.bed"), col_names = FALSE)
  df_de |>
    mutate(name = gene) |>
    filter(ltsr > 0.9, beta < 0) |>
    separate("gene", c("chr", "left", "right"), sep = "-") |>
    dplyr::select(chr, left, right, name, beta) |>
    write_tsv(file = here::here("output", DOCNAME, "HSC_dev_da_lmm.peak_dn.bed"), col_names = FALSE)
}


sessionInfo()
Sys.Date()
