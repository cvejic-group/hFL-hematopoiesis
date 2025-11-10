#!/usr/bin/env Rscript

# nohup Rscript code/01.load_mat_from_bed.R > logs/01.load_mat_from_bed.log &

DOCNAME <- "DNAm"
dir.create(here::here("output", DOCNAME), showWarnings = FALSE)

library(tidyverse)
library(data.table)
library(bsseq)
source(here::here("utils/DNAm.R"))

# donor information
sam_info <- read_csv("/work/DevM_analysis/utils/sam_info.DevM.01.04.25.csv") |>
  dplyr::select(donorID, PCW, Sex) |>
  distinct() |>
  filter(PCW != "Mixed") |>
  mutate(
    PCW = as.numeric(PCW),
    Sex = as.factor(Sex)
  )

# add Nanopore batch
df_batch <- read_tsv(here::here("data/NanoporeBatch.tsv"))
sam_info <- df_batch |>
  left_join(sam_info, by = "donorID") |>
  mutate(SortBatch = factor(SortBatch),
         SampleNumber = factor(SampleNumber)) |>
  column_to_rownames("donorID")
dim(sam_info)
str(sam_info)

# load every sample with bedmethyl.gz available
bed_dir <- "/work/DevM_analysis/data/gDNA/"
samples <- c()
bed_files <- c()
for (d in rownames(sam_info)) {
  f <- paste0(bed_dir, d, "/", d, ".bedmethyl.gz")
  if (file.exists(f)) {
    bed_files <- c(bed_files, f)
    samples <- c(samples, d)
  }
}
length(bed_files)

# Read each bed file with explicit column names
beds <- lapply(bed_files, function(f) {
  dt <- fread(f, header = FALSE)
  setnames(dt, c("chrom", "start", "end", "mod_code",
                 "score", "strand", "compat_start", "compat_end", "color",
                 "N_valid_cov", "frac_mod", "N_mod", "N_canonical",
                 "N_other_mod", "N_delete", "N_fail", "N_diff", "N_nocall"))
  return(dt)
})

# 5mC
bs_5mC  <- build_bsseq_for_mod(sams = samples, bed_lst = beds, mod_code_var = "m")
pData(bs_5mC) <- sam_info[rownames(pData(bs_5mC)),]
bs_5mC

# save
qs2::qs_save(bs_5mC, file = here::here("output", DOCNAME, "bs_5mC_raw.qs2"))

# sam info
xx <- as.data.frame(pData(bs_5mC)) |>
  rownames_to_column("donorID")
saveRDS(xx, file = here::here("output", DOCNAME, "DNAm_samInfo.rds"))
xx |>
  write_tsv(file = here::here("output", DOCNAME, "DNAm_samInfo.tsv"))

# 5hmC
bs_5hmC <- build_bsseq_for_mod(sams = samples, bed_lst = beds, mod_code_var = "h")
pData(bs_5hmC) <- sam_info[rownames(pData(bs_5hmC)),]
bs_5hmC

# save
qs2::qs_save(bs_5hmC, file = here::here("output", DOCNAME, "bs_5hmC_raw.qs2"))

sessionInfo()
Sys.Date()
