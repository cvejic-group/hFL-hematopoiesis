###################################
#                                 #
#    02.1 Joint cell filtering    #
#           ATAC & RNA            #
#                                 #
###################################

## Project: Developmental multiome
## Date: 7.03.2022
## Author: Pernille
## Script name: 02.1.cell_filtering.R

## Script description:
## Cell filtering

##Summary to keep:
## ATAC:           nCount_peaks > 3000
##                 nucleosome_signal < 1.5, TSS.enrichment > 3

## RNA:            nFeatures > 300, nmads < 2 higher Mt, nmads < 2.5 higher rb

.libPaths("/work/home/Software/R/")

library(Seurat)
library(Signac)
library(tidyverse)
library(tidyr)
library(dplyr)

##---------------------------------------------##
##----------------0. Load data-----------------##
##---------------------------------------------##

# sam_info <- read.csv("/work/DevM_analysis/utils/sam_info.DevM.csv", row.names = 1, header = T) %>%
#   filter(GoodRNA == "Yes")
sam_info <- read.csv("/work/DevM_analysis/utils/sam_info.DevM.csv",
                     header = T) %>%
  distinct(Library, .keep_all = TRUE) %>%
  column_to_rownames("Library")
meta_cells <- read.table("/work/DevM_analysis/01.annotation/02.cleandata/data/cell-info.txt", sep = "\t")

#For CB used:
#sam_info_CB <- read.table(file = "/work/DevM_analysis/utils/sam_info.DevM.20.08.24.txt", sep = "\t", header = T) %>%
#  distinct(libraryID, .keep_all = TRUE) %>%
#  column_to_rownames("libraryID") %>%
#  dplyr::filter(Tissue == "CB")

#meta_cells_CB <- read.table("/work/DevM_analysis/01.annotation/02.cleandata/data/CB_cell-info_filtering.txt", sep = "\t")

##---------------------------------------------##
##-----------0. Load RNA doublets--------------##
##---------------------------------------------##
dbl_res <- lapply(rownames(sam_info),
                  FUN = function(s) read.csv(paste0("/work/DevM_analysis/01.annotation/01.doublet/data/", s, ".csv"),
                                             row.names = 1))
names(dbl_res) <- rownames(sam_info)

# Add libraryID and integration_barcode
dbl_res <- lapply(names(dbl_res), FUN = function(s){
  dbl_res[[s]] <- dbl_res[[s]] %>%
    mutate(libraryID = s) %>%
    rownames_to_column("gex_barcode")
  return(dbl_res[[s]])})

dbl_res <- do.call(rbind.data.frame, dbl_res)

meta_cells <- left_join(meta_cells, dbl_res[, c("scDblFinder.class", "gex_barcode", "libraryID")],
                        by = join_by(gex_barcode, libraryID))
meta_cells <- dplyr::rename(meta_cells, "status.scDblFinder" = "scDblFinder.class")

##---------------------------------------------##
##-------0. Load souporcell demultiplex--------##
##---------------------------------------------##

meta_cells$sampleID <- meta_cells$libraryID
meta_cells$donorID <- sapply(strsplit(meta_cells$libraryID, "-"), '[',1)

pool_info <- data.frame(ID = c("FLpool1", "FLpool2", "FLpool3", "FLpool4"),
                        repeated = c(TRUE, FALSE, TRUE, FALSE))


extract_souporcell_results <- function(x){
  sop <- read.table(paste0("/work/DevM_analysis/01.annotation/00.demultiplex/data/", x["ID"], ".tsv"),
                    sep = "\t", header = T)
  if(x["repeated"]){
    sop$libraryID <- sapply(strsplit(sop$barcode, "_"), '[',1)
    sop$barcode <- str_remove_all(sop$barcode, paste0(sop$libraryID, "_"))
  }else{
    sop$libraryID <- x["ID"]
  }

  sop$donorID <- paste0(sapply(strsplit(sop$libraryID, "-"), '[',1), "d", sop$assignment)

  if(!is.na(strsplit(sop$libraryID[1], split = "-")[[1]][2])){
    sop <- sop %>%  mutate("sampleID" = paste0(donorID, "-", str_sub(libraryID, start=-1)))
  }else{
    sop$sampleID <- sop$donorID
  }

  sop <- sop %>% mutate("integration_barcode" = paste(libraryID, barcode, sep = "_")) %>%
    dplyr::select(integration_barcode, status, libraryID, donorID, sampleID)

  return(sop)
}

soc_res <- apply(pool_info, MARGIN = 1, FUN = extract_souporcell_results)
soc_res <- do.call(rbind.data.frame, soc_res)

# Merge with meta
meta_cells <- rquery::natural_join(soc_res[, c("integration_barcode", "status", "sampleID", "donorID", "libraryID")],
                                   meta_cells,
                                   by = "integration_barcode",
                                   jointype = "FULL")
meta_cells <- dplyr::rename(meta_cells, "status.souprorcell" = "status")

##---------------------------------------------##
##-----------1. Define low QC cells------------##
##---------------------------------------------##

meta_cells <- meta_cells %>% mutate(
    HighQualityCell = case_when(
      nCount_RNA < 250 | nFeature_RNA < 300 | nmads2.Mt.higher.log == 1 | nmads2.5.rb.higher.log == 1 |
        nCount_peaks < 3000 | TSS.enrichment < 3 | nucleosome_signal > 1.5 |
        status.scDblFinder == "doublet" |
        status.souprorcell %in% c("doublet", "unassigned") ~ 0,
      TRUE ~ 1
    ))

meta_cells <- meta_cells %>%
  group_by(donorID) %>%
  mutate(nHighQualityCells_perDonorID = sum(HighQualityCell == 1))

meta_cells[meta_cells$nHighQualityCells_perDonorID < 350, "HighQualityCell"] <- 0

write.table(meta_cells, "/work/DevM_analysis/01.annotation/02.cleandata/data/FL_cell-info_filtering.txt", sep = "\t")

##---------------------------------------------##
##------------2. Filtering seu obj-------------##
##---------------------------------------------##

my_data <- list.files(path = "/work/DevM_analysis/01.annotation/02.cleandata/data/seu_obj",
                      full.names = TRUE, recursive = F)

lapply(my_data, FUN = function(path){

  cat(path, "\n")
  seu <- readRDS(path)
  libID <- strsplit(strsplit(path, "/")[[1]][8], "_")[[1]][1]

  seu@meta.data <- seu@meta.data |>
    rownames_to_column("gex_barcode") |>
    right_join(meta_cells_CB %>%
                 dplyr::filter(libraryID == libID) %>%
                 dplyr::select(HighQualityCell, donorID, sampleID, libraryID, integration_barcode, gex_barcode),
               by = "gex_barcode") |>
    column_to_rownames("gex_barcode")

  Idents(seu) <- "HighQualityCell"
  seu <- subset(x = seu, idents = 1)

  if(grepl("rna", path, fixed=TRUE)){
    seu <- seu %>% NormalizeData() %>%
      FindVariableFeatures() %>%
      ScaleData() %>%
      RunPCA() %>%
      RunUMAP(dims = 1:30)
  }
  if(grepl("atac", path, fixed=TRUE)){
    seu <- seu %>% RunTFIDF() %>%
      FindTopFeatures(min.cutoff = 'q25') %>%
      RunSVD(n = 30, reduction.name = "lsi") %>%
      RunUMAP(reduction = "lsi", dims = 2:30)
  }

  saveRDS(seu, path)
})
