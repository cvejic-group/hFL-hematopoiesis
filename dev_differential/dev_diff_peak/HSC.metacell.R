#!/usr/bin/env Rscript

# R432
# cd /work/home/project/20231127_DevM/devm_r432
# nohup Rscript code/HSC.metacell.R > logs/HSC.metacell.log &

DOCNAME <- "HSC.metacell"
dir.create(here::here("output", DOCNAME), showWarnings = FALSE)

library(tidyverse)
library(Seurat)

# metacell group
df_mc <- read_csv("~/project/20231127_DevM/E2G/metacell/hsc_mc_bySample.csv")
length(unique(df_mc$mc_group))
head(df_mc)

# srat
srat <- readRDS("/work/DevM_analysis/01.annotation/11.subclustering/HSC/data/FL_rna_clustered.v01.rds")
srat

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
pm <- readRDS("~/project/20231127_DevM/E2G/SCENT/HSC/pm.rds")
dim(pm)

# nFrags
df_archr <- read_csv("~/project/20231127_DevM/archr_rproj/output/archr_48FL.eachCell/HSC/archr_cellmeta.csv") |>
  mutate(barcode = str_replace(barcode, "#", "_")) |>
  dplyr::select(barcode, nFrags)
# mdata
mdata <- srat@meta.data |>
  rownames_to_column("barcode") |>
  left_join(df_archr, by = "barcode") |>
  left_join(df_mc, by = "barcode") |> # mc_group
  dplyr::select(-celltype, -starts_with("RNA_snn")) |>
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
    m_frag = log(mean(nFrags)),
    libraryID = unique(libraryID),
    donorID = unique(donorID),
    sampleID = unique(sampleID),
    Sex = unique(Sex),
    PCW = unique(PCW),
    Batch = unique(Batch)
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

# log norm
gm_lognorm <- NormalizeData(gm)
pm_lognorm <- NormalizeData(pm)
mc_gm_lognorm <- NormalizeData(mc_gm)
mc_pm_lognorm <- NormalizeData(mc_pm)


# save single-cell
qs::qsave(pm, file = here::here("output", DOCNAME, "pm.qs"))
qs::qsave(gm, file = here::here("output", DOCNAME, "gm.qs"))
qs::qsave(pm_lognorm, file = here::here("output", DOCNAME, "pm_lognorm.qs"))
qs::qsave(gm_lognorm, file = here::here("output", DOCNAME, "gm_lognorm.qs"))
qs::qsave(mdata, file = here::here("output", DOCNAME, "mdata.qs"))

# save metacell
qs::qsave(mc_pm, file = here::here("output", DOCNAME, "mc_pm.qs"))
qs::qsave(mc_gm, file = here::here("output", DOCNAME, "mc_gm.qs"))
qs::qsave(mc_pm_lognorm, file = here::here("output", DOCNAME, "mc_pm_lognorm.qs"))
qs::qsave(mc_gm_lognorm, file = here::here("output", DOCNAME, "mc_gm_lognorm.qs"))
qs::qsave(mc_mdata, file = here::here("output", DOCNAME, "mc_mdata.qs"))


# filtering strategies for genes
# 1 count 1% single-cell (LMM RNA analysis, ~439 cells)
gene_frac_expressed <- rowSums(gm > 0)/ncol(gm)
flag_rna <- gene_frac_expressed > 0.01
qs::qsave(flag_rna, file = here::here("output", DOCNAME, "hsc_gene.flag_1_1.qs"))
qs::qsave(gene_frac_expressed, file = here::here("output", DOCNAME, "hsc_gene_frac_expressed.qs"))


# different filtering strategies for peaks
# 1 count 1% single-cell (LMM RNA analysis, ~439 cells)
peak_frac_expressed <- rowSums(pm > 0)/ncol(pm)
qs::qsave(peak_frac_expressed, file = here::here("output", DOCNAME, "hsc_peak_frac_expressed.qs"))
flag_1_1 <- peak_frac_expressed > 0.01
sum(flag_1_1)
sum(flag_1_1)/nrow(pm)
peak_names <- rownames(pm)[flag_1_1]
data.frame(peak_name = peak_names) |>
  separate(peak_name, c("chrom", "left", "right")) |>
  mutate(peak_name = peak_names) |>
  write_tsv(file = here::here("output", DOCNAME, "hsc_peak.keep_1_1.bed"), col_names = FALSE)

# 1 count 1000 single-cell (~2% cells)
flag_1_1k <- rowSums(pm > 0) > 1000
sum(flag_1_1k)
sum(flag_1_1k)/nrow(pm)
peak_names <- rownames(pm)[flag_1_1k]
data.frame(peak_name = peak_names) |>
  separate(peak_name, c("chrom", "left", "right")) |>
  mutate(peak_name = peak_names) |>
  write_tsv(file = here::here("output", DOCNAME, "hsc_peak.keep_1_1k.bed"), col_names = FALSE)

# 1 count 5% single-cell (SCENT way, ~2196 cells)
flag_1_5 <- peak_frac_expressed > 0.05
sum(flag_1_5)
sum(flag_1_5)/nrow(pm)
peak_names <- rownames(pm)[flag_1_5]
data.frame(peak_name = peak_names) |>
  separate(peak_name, c("chrom", "left", "right")) |>
  mutate(peak_name = peak_names) |>
  write_tsv(file = here::here("output", DOCNAME, "hsc_peak.keep_1_5.bed"), col_names = FALSE)

# save different filtering
qs::qsave(flag_1_1, file = here::here("output", DOCNAME, "hsc_peak.flag_1_1.qs"))
qs::qsave(flag_1_1k, file = here::here("output", DOCNAME, "hsc_peak.flag_1_1k.qs"))
qs::qsave(flag_1_5, file = here::here("output", DOCNAME, "hsc_peak.flag_1_5.qs"))


# scale - mc_pm (1 count 5% single-cell)
srat <- CreateSeuratObject(counts = mc_pm)
srat <- NormalizeData(srat)
srat <- ScaleData(srat, features = rownames(srat)[flag_1_5])
mc_pm_scale = GetAssayData(srat, layer = "scale.data")
qs::qsave(mc_pm_scale, file = here::here("output", DOCNAME, "mc_pm_scale.qs"))

# scale - pm (1 count 5% single-cell)
srat <- CreateSeuratObject(counts = pm)
srat <- NormalizeData(srat)
srat <- ScaleData(srat, features = rownames(srat)[flag_1_5])
pm_scale = GetAssayData(srat, layer = "scale.data")
qs::qsave(pm_scale, file = here::here("output", DOCNAME, "pm_scale.qs"))

# scale - mc_gm (1 count 1% single-cell)
srat <- CreateSeuratObject(counts = mc_gm)
srat <- NormalizeData(srat)
srat <- ScaleData(srat, features = rownames(srat)[flag_rna])
mc_gm_scale = GetAssayData(srat, layer = "scale.data")
qs::qsave(mc_gm_scale, file = here::here("output", DOCNAME, "mc_gm_scale.qs"))

# scale - gm (1 count 1% single-cell)
srat <- CreateSeuratObject(counts = gm)
srat <- NormalizeData(srat)
srat <- ScaleData(srat, features = rownames(srat)[flag_rna])
gm_scale = GetAssayData(srat, layer = "scale.data")
qs::qsave(gm_scale, file = here::here("output", DOCNAME, "gm_scale.qs"))


sessionInfo()
Sys.Date()

# olp with CD34+ HSPC H3K27ac
#27930/149451
#27383/126611
#24125/74167

# olp with scE2G elements
#63804/149451
#61737/126611
#49365/74167

# olp with either
#69798/149451
#67258/126611
#52741/74167


# scE2G elements to here
#45049/50129
#43275/50129
#34447/50129


