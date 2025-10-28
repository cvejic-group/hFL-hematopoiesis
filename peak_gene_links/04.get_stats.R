#!/usr/bin/env Rscript

# R432
# cd ~/project/20231127_DevM/E2G/lmmPG
# nohup Rscript 03.make_bedpe.R > 03.make_bedpe.log &

# set up
.libPaths("~/project/20231127_DevM/devm_r432/renv/library/R-4.3/x86_64-pc-linux-gnu")
library(tidyverse)

quiet <- function(x) suppressMessages(suppressWarnings(x))
get_stats <- function(input_cell=NULL) {
  print(input_cell)
  # set up
  work_dir <- "~/project/20231127_DevM/E2G/lmmPG"
  out_dir <- "/work/DevM_analysis/data/E2G/lmmPG"
  # num gene
  num_gene <- sum(x <- qs::qread(file.path(work_dir, input_cell, "gene.flag_1_1.qs")))
  # metacell num
  num_metacell <- ncol(x <- qs::qread(file.path(work_dir, input_cell, "mc_gm.qs")))
  # num peak
  num_peak <- sum(x <- qs::qread(file.path(work_dir, input_cell, "peak.flag_1_5.qs")))
  # pg link tested
  pg_in <- qs::qread(file.path(work_dir, input_cell, "cis.g2p_list.qs"))
  num_pg_tested <- sum(sapply(pg_in, function(x){nrow(x)}))
  num_sig_pg <- nrow(x <- quiet(read_tsv(file.path(out_dir, paste0(input_cell, ".lmmPG_FDR20.bedpe")), col_names = F)))
  return(
    data.frame(
      celltype = input_cell,
      num_metacell = num_metacell,
      num_gene = num_gene,
      num_peak = num_peak,
      num_pg_tested = num_pg_tested,
      num_pg_sig = num_sig_pg
    )
  )
}

# run
cells <- c("GP", "Granulocyte",
           "MEMP-t", "MEMP", "MEP", "MEMP-Mast-Ery", "MEMP-Ery", "Early-Ery", "Late-Ery",
           "MEMP-MK", "MK", "MastP-t", "MastP", "Mast",
           "MDP", "Monocyte", "Kupffer", "cDC1", "cDC2", "pDC", "ASDC",
           "LMPP", "LP", "Cycling-LP", "PreProB", "ProB-1", "ProB-2", "Large-PreB", "Small-PreB", "IM-B",
           "NK", "ILCP", "T")
tmp_lst <- lapply(cells, function(i){get_stats(input_cell = i)})
df <- do.call(rbind, tmp_lst)

# num_cell
cell_num <- read_csv("/work/DevM_analysis/01.annotation/11.subclustering/blood/data/FL_wnn_cellmeta.v00.csv") |>
  group_by(anno_wnn_v51) |>
  summarise(num_cell = n()) |>
  dplyr::rename(celltype = anno_wnn_v51)
# out
df <- df |>
  left_join(cell_num, by = "celltype")
df |>
  write_csv(file = "pg_stats.csv")

