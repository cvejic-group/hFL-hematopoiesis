###################################
#                                 #
#      02.0 Prepare objects       # 
#          & QC metrics           #
#                                 #
###################################

## Project: Developmental multiome
## Author: Pernille 
## Script name: 02.0.prepare_data.R

## Script description: 
## Calculate QC metrics and create individual Seurat objects for each library 

.libPaths("/work/home/Software/R/")

library(Seurat)#; library(Signac)
library(scater); library(readxl)
library(tidyverse)
library(dplyr)
library(EnsDb.Hsapiens.v86); library(GenomeInfoDb)


cr_path <- "/work/cellranger-arc/"
meta_cells <- list()

##---------------------------------------------##
##---------------0. Loading data---------------##
##---------------------------------------------##

## ------- 0. sample info
# FL samples were run from: 
sam_info <- read.table(file = "/work/DevM_analysis/utils/sam_info.DevM.12.08.24.txt", sep = "\t", header = T) %>% 
  distinct(libraryID, .keep_all = TRUE) 
rownames(sam_info) <- sam_info$libraryID

#CB samples were 
sam_info_CB <- read.table(file = "/work/DevM_analysis/utils/sam_info.DevM.20.08.24.txt", sep = "\t", header = T) %>% 
  distinct(libraryID, .keep_all = TRUE) 
rownames(sam_info_CB) <- sam_info_CB$libraryID
sam_info_CB <- sam_info_CB %>% dplyr::filter(Tissue == "CB")

## ------- 0. Annotation genome
library(AnnotationHub)
ah <- AnnotationHub()
query(AnnotationHub(), c("Homo", "Ensdb", "98"))
ensdb <- ah[["AH75011"]]
annotations <- GetGRangesFromEnsDb(ensdb = ensdb,
                                   standard.chromosomes = T)
seqlevelsStyle(annotations) <- "UCSC"

## ------- 1. Read data (per library)
#for(s in rownames(sam_info)){
for(s in rownames(sam_info_CB)){ 
  cat(s, "starts", "\n")
  
  input <- Read10X_h5(paste0(cr_path, sam_info_CB[ s, "CellRanger"], "/filtered_feature_bc_matrix.h5"))
  
  rna <- input$`Gene Expression`
  atac <- input$Peaks
  
  rm(input)
  
  ##---------------------------------------------##
  ##----------------1. RNA-seq-------------------##
  ##---------------------------------------------##
  
  cat("  rna starts", "\n")
  
  ## ------- 0. create seu obj
  rna <- CreateSeuratObject(rna)
  
  ## ------- 1. Calculate QC metrics
  RPS.RPL.genes <- c(rownames(rna)[grep("^RPS", rownames(rna))],
                     rownames(rna)[grep("^RPL", rownames(rna))])
  rna$percent.mt <- PercentageFeatureSet(rna, pattern = "^MT-")
  rna$percent.rb <- PercentageFeatureSet(rna, features = RPS.RPL.genes)
  
  ## ------- 2. Calculate QC metrics ignoring Mt genes
  rna$nCount_excludingMt <- colSums(rna[["RNA"]]["counts"][-c(grep("^MT-", rownames(rna))),])
  rna$nFeature_excludingMt <- colSums(rna[["RNA"]]["counts"][-c(grep("^MT-", rownames(rna))),] > 0)
  
  rna$nCount_Mt <- colSums(GetAssayData(rna, slot = "counts")[c(grep("^MT-", rownames(rna))),])
  
  rna$nmads2.Mt.higher.log <- scater::isOutlier(rna$percent.mt, nmads = 2, type = "higher", log = T)
  rna$nmads2.5.rb.higher.log <- scater::isOutlier(rna$percent.rb, nmads = 2.5, type = "higher", log = T)
  
  ##---------------------------------------------##
  ##----------------2. ATAC-seq------------------##
  ##---------------------------------------------##
  
  cat("  atac starts", "\n")
  
  ## ------- 0. Create chrom assay, seu obj, add meta
  atac <- CreateChromatinAssay(
    counts = atac,
    sep = c(":", "-"),
    genome = 'hg38',
    fragments = paste0(cr_path, sam_info_CB[ s, "CellRanger"], "/atac_fragments.tsv.gz"),
    annotation = annotations,
    min.cells = 10
  )
  atac <- CreateSeuratObject(counts = atac,
                             assay = "peaks")
  
  ## ------- 1. Calculate QC metrics
  atac <- NucleosomeSignal(object = atac)
  atac <- TSSEnrichment(atac, fast = T)
  
  meta_cr <- read.csv(file = paste0(cr_path, sam_info_CB[ s, "CellRanger"], "/per_barcode_metrics.csv"),
                      header = TRUE,row.names = 1)
  atac@meta.data = atac@meta.data |>
    rownames_to_column("gex_barcode") |>
    left_join(y = meta_cr, by = "gex_barcode") |>
    column_to_rownames("gex_barcode")
  
  atac$pct_reads_in_peaks <- atac$atac_peak_region_fragments / atac$atac_fragments * 100
  atac$pct_mt <- atac$atac_mitochondrial_reads/atac$atac_raw_reads *100
  
  ##---------------------------------------------##
  ##------------------3. meta--------------------##
  ##---------------------------------------------##
  meta_cells[[s]] <- left_join(x = rna@meta.data %>% 
                                 rownames_to_column("gex_barcode"),
                               y = atac@meta.data %>% 
                                 rownames_to_column("gex_barcode"), 
                               by = "gex_barcode")
  meta_cells[[s]]$libraryID <- s
  
  ##---------------------------------------------##
  ##------------------4. Save--------------------##
  ##---------------------------------------------##
  cat(" ", s,  "is being saved", "\n")
  
  saveRDS(object = rna,
          paste0("/work/DevM_analysis/01.annotation/02.cleandata/data/seu_obj/", s, "_rna_seu.rds"))
  saveRDS(object = atac,
          paste0("/work/DevM_analysis/01.annotation/02.cleandata/data/seu_obj/", s, "_atac_seu.rds"))
  
  cat(" done", "\n")
}

meta_cells <- do.call(rbind.data.frame, meta_cells)
meta_cells$integration_barcode <- paste0(meta_cells$libraryID, "_", meta_cells$gex_barcode)
rownames(meta_cells) <- meta_cells$integration_barcode

meta_cells_FL <- left_join(meta_cells, sam_info %>% select("PCW", "libraryID"), by = "libraryID")
meta_cells_CB <- left_join(meta_cells, sam_info_CB[, c("PCW", "libraryID")], by = "libraryID")

write.table(meta_cells_FL, "/work/DevM_analysis/01.annotation/02.cleandata/data/FL_cell-info_filtering.txt", sep = "\t")
write.table(meta_cells_CB, "/work/DevM_analysis/01.annotation/02.cleandata/data/CB_cell-info_filtering.txt", sep = "\t")
