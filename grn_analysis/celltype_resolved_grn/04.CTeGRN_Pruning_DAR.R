#######################
### Pruning CT-eGRN ###
#######################

###-------------------Note-------------------###
###-------      Run with R-4.4.2      -------###
###------------------------------------------###

##################
# Initialization #
##################

setwd("~/local_data/proj/Dev_Multiome/04.regulome_R/06.SCENICplus_CTeGRN/")
source("./00.Initialization.R")

library(AUCell)
library(ArchR)
library(caret)
library(pROC)
library(data.table)
library(future)
library(furrr)


##########################
### Load Filtered eGRN ###
##########################

# Load SCENIC+ eRegulon metadata
raw_eRegulon_meta.df <- readRDS(paste0(DATA_DIR, "eRegulon_raw_meta.rds"))
BA_eRegulon_meta.df <- raw_eRegulon_meta.df %>% dplyr::filter(class == "+/+")

# Form CT-eRegulons list
CT_eRegulon.l <-readRDS(paste0(RES_DIR, 
                               "eRegulons_CT_Filter/CT_eGRN_Filter_Metadata.rds"))

################################ 
### Pruning CT-eGRN with DAR ###
################################

# Import DARs
DAR.path <- "~/local_data/proj/Dev_Multiome/04.regulome/scp_ALL_pcw/lumi_outs/DAR_outs/scp_snakemake_ALL_pcw/data/region_sets/DARs_cell_type/"
DAR.peaks.l <- list()
for (CT in CT_ORDER) {
  DAR.peaks.l[[CT]] <- rtracklayer::import(paste0(DAR.path, CT, ".bed"))
}

# Prune eRegulon according to DARs
### Function
CT_eRegulon_prune <- function(CT_eRegulon_filter_meta, DAR.peaks.gr, CT.idx){
  CT_eRegulon_prune_meta <- data.frame()
  discard.idx <- 0
  for (TF.idx in unique(CT_eRegulon_filter_meta$TF)) {
    eRegulon.ori <- CT_eRegulon_filter_meta %>% filter(TF == TF.idx)
    TR.gr <- as(eRegulon.ori$Region, "GRanges")
    TR.overlap <- findOverlaps(TR.gr, DAR.peaks.gr)
    if (is_empty(TR.overlap)) {
      eRegulon.prune <- NA
      discard.idx <- discard.idx + 1
    }else{
      eRegulon.prune <- eRegulon.ori[TR.overlap@from ,]
      eRegulon.prune$eRegulon_CT <- paste0(eRegulon.prune$TF,
                                           "_", CT.idx,
                                           "_", eRegulon.prune$class)
      eRegulon.prune$TG_num <- length(unique(eRegulon.prune$Gene))
      eRegulon.prune$Region_num <- length(unique(eRegulon.prune$Region))
      if (unique(eRegulon.prune$TG_num) < 0.1 * unique(eRegulon.ori$TG_num) | unique(eRegulon.prune$TG_num) < 5) {
        eRegulon.prune <- NA
        discard.idx <- discard.idx + 1
      }else{
        CT_eRegulon_prune_meta <- rbind(CT_eRegulon_prune_meta, eRegulon.prune)
      }
    }
  }
  return(list(CT_eRegulon_prune_meta = CT_eRegulon_prune_meta,
              discard.rate = discard.idx/length(unique(CT_eRegulon_filter_meta$TF)),
              eRegulon_ori = length(unique(CT_eRegulon_filter_meta$TF)),
              eRegulon_purne = length(unique(CT_eRegulon_prune_meta$TF))))
}

# Prune eRegulons
prune_meta.l <- list()
for (CT in CT_ORDER) {
  DAR.peaks.gr <- DAR.peaks.l[[CT]]
  CT_eRegulon_filter_meta <- CT_eRegulon.l[[CT]]
  prune_meta.l[[CT]] <- CT_eRegulon_prune(CT_eRegulon_filter_meta, DAR.peaks.gr, CT)
}
# Save Pruned results
saveRDS(prune_meta.l, 
        paste0(RES_DIR, "eRegulons_CT_Prune/CT_eGRN_DAR_Prune_Metadata.rds"))
### Save to xlsx
library(openxlsx)
wb <- createWorkbook()
for(CT in CT_ORDER){
  addWorksheet(wb, CT)
  writeData(wb, CT, prune_meta.l[[CT]]$CT_eRegulon_prune_meta)
}
saveWorkbook(wb, paste0(RES_DIR, "eRegulons_CT_Prune/CT_eGRN_DAR_Prune_Metadata.xlsx"), 
             overwrite = TRUE)
