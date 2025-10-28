#!/usr/bin/env Rscript

# R432
# cd ~/project/20231127_DevM/E2G/lmmPG
# nohup Rscript 03.make_bedpe.R > 03.make_bedpe.log &

# set up
.libPaths("~/project/20231127_DevM/devm_r432/renv/library/R-4.3/x86_64-pc-linux-gnu")
library(tidyverse)

quiet <- function(x) suppressMessages(suppressWarnings(x))
make_bedpe <- function(input_cell=NULL, gene_loc=NULL) {
  print(input_cell)
  # set up
  work_dir <- "~/project/20231127_DevM/E2G/lmmPG"
  out_dir <- "/work/DevM_analysis/data/E2G/lmmPG"
  # load
  fname <- file.path(work_dir, input_cell, paste0(input_cell, ".lmmPG.tsv"))
  df_raw <- quiet(read_delim(fname, delim = " ")) |>
    mutate(padj = p.adjust(p_value, method = "BH"),
           EG = paste(peak, gene, sep = ":")) |>
    filter(beta_atac > 0)
  # bedpe
  out_file <- file.path(out_dir, paste0(input_cell, ".lmmPG_FDR20.bedpe"))
  df <- df_raw |> filter(padj < 0.2) # same as HSC
  print(nrow(df))
  df |>
    left_join(gene_loc, by = "gene") |>
    separate(col = "peak", into = c("chr1", "left1", "right1"), sep = "-") |>
    filter(chr1 == chr2) |>
    dplyr::select(chr1, left1, right1, chr2, left2, right2, EG, beta_atac) |>
    mutate(strand1 = ".", strand2 = ".") |>
    mutate(chr1 = factor(chr1, levels = paste0("chr", c(1:22, "X")))) |>
    arrange(chr1, left1) |>
    write_tsv(
      file = out_file,
      col_names = FALSE
    )
}

# gene loc
genebed_file <- "~/RefData/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/GeneBody_500kb_margin.bed"
gene_loc <- quiet(read_tsv(genebed_file, col_names = FALSE)) |>
  dplyr::select(chr2 = X1, left2 = X2, right2 = X3, gene = X4)

# test
#i <- "GP"
#make_bedpe(input_cell = i, gene_loc = gene_loc)

# run
cells <- c("GP", "Granulocyte",
           "MEMP-t", "MEMP", "MEP", "MEMP-Mast-Ery", "MEMP-Ery", "Early-Ery", "Late-Ery",
           "MEMP-MK", "MK", "MastP-t", "MastP", "Mast",
           "MDP", "Monocyte", "Kupffer", "cDC1", "cDC2", "pDC", "ASDC",
           "LMPP", "LP", "Cycling-LP", "PreProB", "ProB-1", "ProB-2", "Large-PreB", "Small-PreB", "IM-B",
           "NK", "ILCP", "T")
for (i in cells) {
  make_bedpe(input_cell = i, gene_loc = gene_loc)
}

devtools::session_info()
Sys.Date()
