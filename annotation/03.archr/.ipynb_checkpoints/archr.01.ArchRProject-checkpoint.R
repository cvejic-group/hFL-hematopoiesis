#!/usr/bin/env Rscript

# run
# nohup Rscript analysis/archr.01.ArchRProject.R > analysis/archr.01.ArchRProject.log &

DOCNAME <- "archr_48FL"
dir.create(here::here("output", DOCNAME), showWarnings = FALSE)

library(tidyverse)
library(ArchR)
library(Seurat)


# set up
set.seed(1)
# default number of threads
addArchRThreads(threads = 12)
# add a reference genome annotation for ArchR
addArchRGenome("hg38")
# chr prefix to exclude KI/GL scaffolds
addArchRChrPrefix(chrPrefix = TRUE)
# cor cutoff
corCutOff = 0.5
# chr to exclude
excludeChr = c("chrM")

# load samples
sam_info <- read_tsv(here::here("data/sam_info.tsv")) |>
  dplyr::select(libraryID, CellRanger) |>
  distinct()
sample_names <- sam_info$libraryID

# cellranger-arc dir
cr_dir <- "/work/home/data/cellranger-arc/"
inpu_files <- paste0(cr_dir, sam_info$CellRanger, "/atac_fragments.tsv.gz")

# clean barcodes
bc_dir <- "/work/DevM_analysis/01.annotation/10.integration_joint_clean/data/clean_barcodes_per_lib/"
clean_barcodes <- paste0(bc_dir, sample_names, '.csv')

# out prefix
dir.create(here::here("output", "archr.ArrowFiles"), showWarnings = FALSE)
output_names <- paste0(here::here("output", "archr.ArrowFiles/"), sample_names)

# create
ArrowFiles <- createArrowFiles(
  inputFiles = inpu_files,
  sampleNames = sample_names,
  validBarcodes = getValidBarcodes(clean_barcodes, sample_names),
  outputNames = output_names,
  QCDir = here::here("output", "archr.QualityControl"),
  excludeChr = excludeChr,
  minTSS = 0,
  minFrags = 0,
  maxFrags = Inf,
  minFragSize = 0,
  maxFragSize = Inf,
  addTileMat = TRUE,
  TileMatParams = list(excludeChr = excludeChr),
  addGeneScoreMat = TRUE,
  GeneScoreMatParams = list(excludeChr = excludeChr),
  offsetPlus = 0, offsetMinus = 0,
  force = FALSE
)
ArrowFiles

# Creating an ArchRProject
proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = here::here("output", DOCNAME),
  copyArrows = TRUE # This is recommened so that you maintain an unaltered copy for later usage.
)

# add gene exp.
srat <- readRDS('/work/DevM_analysis/01.annotation/10.integration_joint_clean/data/FL_rna_clustered.v00.rds')
srat <- DietSeurat(srat, counts = TRUE)
# change cell name: _ to #
srat <- RenameCells(srat, new.names = str_replace(colnames(srat), "_", "#"))
seRNA <- as.SingleCellExperiment(srat)
# add rowrange
dummy_seRNA <- import10xFeatureMatrix(
  input = "/work/home/data/cellranger-arc/cellranger-arc202_count_FLpool1_GRCh38-2020-A-2_0_0/filtered_feature_bc_matrix.h5",
  names = c("dummy")
)
rowRanges(seRNA) <- rowRanges(dummy_seRNA)
proj <- addGeneExpressionMatrix(
  input = proj,
  seRNA = seRNA,
  excludeChr = excludeChr
)
# save this one for SCARlink
saveRDS(srat, file = here::here("output", DOCNAME, "scrna_seurat.rds"))

# check
getAvailableMatrices(proj)

# add Donor/Sex/PCW, WNN v05 anno
cellmeta <- read.csv("/work/DevM_analysis/01.annotation/10.integration_joint_clean/data/FL_wnn_cellmeta.v00.csv") |>
  dplyr::select(barcode=X, celltype = anno_wnn_v5, libraryID, donorID, sampleID, PCW, Sex = rna.Sex) |>
  mutate(barcode = str_replace(barcode, "_", "#")) |>
  column_to_rownames("barcode")
cellmeta <- cellmeta[proj$cellNames,]
proj$libraryID <- cellmeta$libraryID
proj$donorID <- cellmeta$donorID
proj$sampleID <- cellmeta$sampleID
proj$PCW <- cellmeta$PCW
proj$Sex <- cellmeta$Sex
proj$celltype <- cellmeta$celltype
proj$celltype <- factor(proj$celltype, levels = c("HSC", "35-MDP?", "GP", "Granulocyte",
                                                  "49-MEMP", "MEMP", "MEP", "MEMP-Mast-Ery", "MEMP-Ery", "Early-Ery", "Late-Ery",
                                                  "MEMP-MK", "MK", "49:2-MastP", "MastP", "Mast",
                                                  "Monocyte", "Kupffer", "cDC1", "cDC2", "pDC", "ASDC",
                                                  "LMPP", "LP", "Cycling_LP", "PreProB", "ProB-1", "ProB-2",
                                                  "Large-PreB", "Small-PreB", "IM-B",
                                                  "NK", "ILCP", "T",
                                                  "Hepatocyte", "Endothelia"))

# https://github.com/GreenleafLab/ArchR/issues/1185
proj$celltype <- str_replace(proj$celltype, '[?]', "")

# save
saveArchRProject(ArchRProj = proj, load = FALSE)

# not binarized tile mat for SCARlink
saveArchRProject(ArchRProj = proj, load = FALSE, outputDirectory = here::here("output", "archr_48FL.nonBinarizedTileMat"))
proj2 <-loadArchRProject(path = here::here('output', "archr_48FL.nonBinarizedTileMat"))
proj2 <- addTileMatrix(input=proj2, binarize=FALSE, tileSize = 500, force=TRUE)
saveArchRProject(ArchRProj = proj2, load = FALSE)

Sys.Date()
sessionInfo()


