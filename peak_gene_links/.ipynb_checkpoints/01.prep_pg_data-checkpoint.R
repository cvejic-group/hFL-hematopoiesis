#!/usr/bin/env Rscript

# R432
# cd ~/project/20231127_DevM/E2G/glmPG
# nohup Rscript 01.prep_pg_data.R > 01.prep_pg_data.log &

# set up
.libPaths("~/project/20231127_DevM/devm_r432/renv/library/R-4.3/x86_64-pc-linux-gnu")
library(tidyverse)
library(Seurat)

prep_mc <- function(input_srat=NULL, input_cell=NULL, df_mc=NULL, df_archr=NULL){
  # set up
  work_dir = "~/project/20231127_DevM/E2G/glmPG"
  out_dir = file.path(work_dir, input_cell)
  dir.create(path = out_dir, showWarnings = FALSE, recursive = TRUE)

  srat <- subset(input_srat, subset = anno_wnn_v51 == input_cell)
  print(srat)

  # no PCW == "Mixed" !!!
  # no PCW == "Mixed" !!!
  # no PCW == "Mixed" !!!
  srat <- subset(srat, subset = PCW != "Mixed")
  srat@meta.data <- droplevels(srat@meta.data)
  table(srat$PCW, useNA = "ifany")

  # GM
  gm <- srat@assays$RNA$counts
  dim(gm)

  # PM
  pm <- readRDS(file.path("~/project/20231127_DevM/E2G/SCENT", input_cell, "pm.rds"))
  dim(pm)

  # mdata
  mdata <- srat@meta.data |>
    rownames_to_column("barcode") |>
    left_join(df_archr, by = "barcode") |>
    left_join(df_mc, by = "barcode") |> # mc_group
    mutate(
      donorID = factor(donorID),
      Sex = factor(Sex),
      Batch = factor(Batch),
      Phase = factor(Phase, levels = c("G1", "S", "G2M"))
    ) |>
    droplevels() |>
    column_to_rownames("barcode") |>
    as.data.frame()
  dim(mdata)
  table(mdata$PCW, useNA = "ifany")
  str(mdata)

  # neat
  pm <- pm[,rownames(mdata)]
  dim(pm)

  # check
  stopifnot(all(colnames(gm) == rownames(mdata)))
  stopifnot(all(colnames(pm) == rownames(mdata)))

  # aggr mdata by mc_group
  mc_mdata <- mdata |>
    rownames_to_column("barcode") |>
    group_by(mc_group) |>
    reframe(
      m_umi = log(mean(nCount_RNA)),
      m_gene = log(mean(nFeature_RNA)),
      m_mt = mean(percent.mt),
      m_rb = mean(percent.rb),
      m_frag = log(mean(nFrags))
    ) |>
    droplevels() |>
    as.data.frame()
  dim(mc_mdata)
  table(mc_mdata$PCW, useNA = "ifany")
  str(mc_mdata)
  # aggr counts
  cut_mat <- sapply(unique(mdata$mc_group), function(g) as.numeric(mdata$mc_group == g))
  cut_mat <- Matrix::Matrix(cut_mat, sparse = TRUE)
  # Summarize by sum
  mc_gm <- gm %*% cut_mat
  mc_pm <- pm %*% cut_mat
  dim(mc_gm)
  dim(mc_pm)
  # change order
  mc_gm = mc_gm[,mc_mdata$mc_group]
  mc_pm = mc_pm[,mc_mdata$mc_group]

  # save mdata
  qs::qsave(mc_mdata, file = file.path(out_dir, "mc_mdata.qs"))
  qs::qsave(mc_gm, file = file.path(out_dir, "mc_gm.qs"))
  qs::qsave(mc_pm, file = file.path(out_dir, "mc_pm.qs"))

  # filtering genes
  # 1 count 1% single-cell
  gene_frac_expressed <- rowSums(gm > 0)/ncol(gm)
  flag_rna <- gene_frac_expressed > 0.01
  qs::qsave(flag_rna, file = file.path(out_dir, "gene.flag_1_1.qs"))
  qs::qsave(gene_frac_expressed, file = file.path(out_dir, "gene_frac_expressed.qs"))

  # filtering peaks
  # 1 count 5% single-cell (SCENT way)
  peak_frac_expressed <- rowSums(pm > 0)/ncol(pm)
  flag_1_5 <- peak_frac_expressed > 0.05
  sum(flag_1_5)
  sum(flag_1_5)/nrow(pm)
  qs::qsave(flag_1_5, file = file.path(out_dir, "peak.flag_1_5.qs"))

  # scale - mc_pm (1 count 5% single-cell)
  srat <- CreateSeuratObject(counts = mc_pm)
  srat <- NormalizeData(srat, verbose = FALSE)
  srat <- ScaleData(srat, features = rownames(srat)[flag_1_5], verbose = FALSE)
  mc_pm_scale <- GetAssayData(srat, layer = "scale.data")
  qs::qsave(mc_pm_scale, file = file.path(out_dir, "mc_pm_scale.qs"))

  # scale - mc_gm (1 count 1% single-cell)
  srat <- CreateSeuratObject(counts = mc_gm)
  srat <- NormalizeData(srat, verbose = FALSE)
  srat <- ScaleData(srat, features = rownames(srat)[flag_rna], verbose = FALSE)
  mc_gm_scale <- GetAssayData(srat, layer = "scale.data")
  qs::qsave(mc_gm_scale, file = file.path(out_dir, "mc_gm_scale.qs"))

  # return
  return(list(mdata = mc_mdata, gm = mc_gm_scale, pm = mc_pm_scale))
}

CreatePeakToGeneList <- function(gm = NULL, pm = NULL,
                                 nbatch = NULL,
                                 genebed = "~/RefData/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/GeneBody_500kb_margin.bed",
                                 tmpfile = "./temporary_atac_peak.bed",
                                 intersectedfile = "./temporary_atac_peak_intersected.bed.gz") {
  peaknames <- rownames(pm) # peak by cell matrix
  peak_bed <- data.frame(
    chr = str_split_fixed(peaknames, "-", 3)[, 1],
    start = str_split_fixed(peaknames, "-", 3)[, 2],
    end = str_split_fixed(peaknames, "-", 3)[, 3],
    peak = peaknames
  )
  write.table(peak_bed, tmpfile, quote = F, row = F, col = F, sep = "\t")
  system(paste("~/bin/bedtools intersect -a", genebed, "-b ", tmpfile, " -wa -wb -loj | gzip -c >", intersectedfile))
  system(paste("rm ", tmpfile)) # delete in case it's used by others
  d <- data.table::fread(intersectedfile, sep = "\t")
  system(paste("rm ", intersectedfile)) # delete in case it's used by others
  d <- data.frame(d)
  d <- d[d$V5 != ".", ]

  # Obtain gene to peak pairs.
  cis.g2p <- d[c("V4", "V8")]
  colnames(cis.g2p) <- c("gene", "peak")
  genes_in_rna <- rownames(gm) # gene by cell matrix
  cis.g2p <- cis.g2p[cis.g2p$gene %in% genes_in_rna, ] # make sure g2p genes are all included in rna matrix
  print(dim(cis.g2p))

  # into small pieces
  cis.g2p$index <- 1:nrow(cis.g2p)
  cis.g2p$batch_index <- Hmisc::cut2(cis.g2p$index, g = nbatch, levels.mean = TRUE)
  cis.g2p_list <- split(cis.g2p, f = cis.g2p$batch_index)
  cis.g2p_list <- lapply(cis.g2p_list, function(x) x[(names(x) %in% c("peak", "gene"))])
  names(cis.g2p_list) <- 1:length(cis.g2p_list)
  return(cis.g2p_list)
}

prep_pg <- function(input_cell=NULL, mdata=NULL, gm=NULL, pm=NULL){

  # set up
  work_dir <- "~/project/20231127_DevM/E2G/glmPG"
  out_dir <- file.path(work_dir, input_cell)
  dir.create(path = out_dir, showWarnings = FALSE, recursive = TRUE)

  stopifnot(all(colnames(gm) == mdata$mc_group))
  stopifnot(all(colnames(pm) == mdata$mc_group))
  # create peak-gene links
  # nbatch to 100 to comply with LUMI limitation about jobs
  cis.g2p_list <- CreatePeakToGeneList(gm = gm, pm = pm, nbatch = 100)
  # save
  qs::qsave(cis.g2p_list, file = file.path(out_dir, "cis.g2p_list.qs"))

  # save list
  obj_lst <- list(
    mdata = mdata,
    rna = gm,
    atac = pm,
    link_lst = cis.g2p_list
  )
  qs::qsave(obj_lst, file = file.path(out_dir, "obj_lst.qs"))
}

# srat
srat_blood <- readRDS("/work/DevM_analysis/01.annotation/11.subclustering/blood/data/FL_rna_clustered.v00.rds")

# metacell group
df_mc <- read_csv("~/project/20231127_DevM/E2G/metacell/blood_celltype_mc.csv")
head(df_mc)

# nFrags
df_archr <- read_csv("~/project/20231127_DevM/archr_rproj/output/archr_48FL/archr_cellmeta.csv") |>
  mutate(barcode = str_replace(barcode, "#", "_")) |>
  dplyr::select(barcode, nFrags)

# loop
for (i in levels(srat_blood$anno_wnn_v51)) {
  print(i)
  tmp_lst <- prep_mc(input_srat=srat_blood, input_cell=i, df_mc=df_mc, df_archr=df_archr)
  x <- prep_pg(input_cell=i, mdata=tmp_lst$mdata, gm=tmp_lst$gm, pm=tmp_lst$pm)
}


sessionInfo()
Sys.Date()
