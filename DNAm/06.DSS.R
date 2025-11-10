#!/usr/bin/env Rscript

# nohup Rscript code/DNAm.DSS.R > logs/DNAm.DSS.log &

DOCNAME <- "DNAm.DSS"
dir.create(here::here("output", DOCNAME), showWarnings = FALSE)

library(tidyverse)
library(data.table)
library(bsseq)
library(DSS)

# load autosomes
bs_5mC_filt <- qs2::qs_read(here::here("output", "DNAm", "bs_5mC_filt.autosomes.qs2"))
bs_5mC_filt

# drop levels
pData(bs_5mC_filt) <- droplevels(pData(bs_5mC_filt))

# fit
DMLfit <- DMLfit.multiFactor(
  bs_5mC_filt,
  design = data.frame(
    PCW = pData(bs_5mC_filt)$PCW,
    Sex  = pData(bs_5mC_filt)$Sex,
    Batch = pData(bs_5mC_filt)$SortBatch
  ),
  formula = ~ PCW + Sex + Batch
)
colnames(DMLfit$X)

# test
DMLtest <- DMLtest.multiFactor(DMLfit, coef = "PCW")
ix <- sort(DMLtest[,"pvals"], index.return=TRUE)$ix
#DMLtest <- DMLtest[ix,]
head(DMLtest[ix,])

# save
qs2::qs_save(DMLfit, file = here::here("output", DOCNAME, "bs_5mC.DSS_DMLfit.qs2"))
qs2::qs_save(DMLtest, file = here::here("output", DOCNAME, "bs_5mC.DSS_DMLtest.qs2"))

# dml to bed
as.data.frame(DMLtest) |>
  mutate(start = pos - 1,
         name = paste(chr, pos, sep = "-")) |>
  dplyr::select(chr, start, end = pos, name, stat, pvals) |>
  write_tsv(file = here::here("output", DOCNAME, "bs_5mC.DSS_DMLtest.bed"), col_names = FALSE)

# dml to bed, fdrs < 0.05
as.data.frame(DMLtest) |>
  filter(fdrs < 0.05) |>
  mutate(start = pos - 1,
         name = paste(chr, pos, sep = "-")) |>
  dplyr::select(chr, start, end = pos, name, stat, fdrs) |>
  write_tsv(file = here::here("output", DOCNAME, "bs_5mC.DSS_DML_fdr005.bed"), col_names = FALSE)
as.data.frame(DMLtest) |>
  filter(fdrs < 0.05, stat > 0) |>
  mutate(start = pos - 1,
         name = paste(chr, pos, sep = "-")) |>
  dplyr::select(chr, start, end = pos, name, stat, fdrs) |>
  write_tsv(file = here::here("output", DOCNAME, "bs_5mC.DSS_DML_fdr005_Up.bed"), col_names = FALSE)
as.data.frame(DMLtest) |>
  filter(fdrs < 0.05, stat < 0) |>
  mutate(start = pos - 1,
         name = paste(chr, pos, sep = "-")) |>
  dplyr::select(chr, start, end = pos, name, stat, fdrs) |>
  write_tsv(file = here::here("output", DOCNAME, "bs_5mC.DSS_DML_fdr005_Dn.bed"), col_names = FALSE)


# dmr
dmrs <- callDMR(DMLtest, p.threshold = 0.05)
dmrs <- dmrs |>
  dplyr::arrange(-areaStat) |>
  mutate(name = paste(chr, start, end, sep = "-")) |>
  dplyr::select(chr, start, end, name, length, nCG, areaStat)

# save
qs2::qs_save(dmrs, file = here::here("output", DOCNAME, "bs_5mC.DSS_DMR_p005.qs2"))

# dmrs to bed
dmrs |>
  write_tsv(file = here::here("output", DOCNAME, "bs_5mC.DSS_DMR_p005.bed"), col_names = FALSE)
dmrs |>
  filter(areaStat > 0) |>
  write_tsv(file = here::here("output", DOCNAME, "bs_5mC.DSS_DMR_p005_Up.bed"), col_names = FALSE)
dmrs |>
  filter(areaStat < 0) |>
  write_tsv(file = here::here("output", DOCNAME, "bs_5mC.DSS_DMR_p005_Dn.bed"), col_names = FALSE)

sessionInfo()
Sys.Date()


