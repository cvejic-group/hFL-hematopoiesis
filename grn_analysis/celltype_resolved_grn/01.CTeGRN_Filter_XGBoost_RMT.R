################################################################
### Using XGBoost RMT to filter Cell-type specific eRegulons ###
################################################################

###-------------------Note-------------------###
###-------      Run with R-4.4.2      -------###
###------------------------------------------###

##################
# Initialization #
##################

setwd("~/local_data/proj/Dev_Multiome/04.regulome_R/01.SCENICplus/04.SCENICplus_CTeGRN/")
source("./00.Initialization.R")

# Load additional packages
## Install packages if not
# packages <- c("xgboost", "caret", "SHAPforxgboost", "dplyr", "pROC", "ggplot2", "data.table")
# install_if_missing <- function(pkg) if (!require(pkg, character.only = TRUE)) BiocManager::install(pkg, lib=.libPaths()[1])
# lapply(packages, install_if_missing)

library(AUCell)
library(ArchR)
library(xgboost)
library(caret)
library(pROC)
library(data.table)
library(SHAPforxgboost)
library(future)
library(furrr)

#################
### Load data ###
#################

# snRNA-seq data
FL.SeuratObj <- readRDS("~/local_data/proj/Dev_Multiome/data/FL_scrna_seurat_20251014.rds")
## Cell Metadata
cell_metadata <- FL.SeuratObj@meta.data

# eRegulon metadata to select BA eRegulons
raw_eRegulon_meta.df <- readRDS(paste0(DATA_DIR, "eRegulon_raw_meta.rds"))
BA_eRegulon_meta.df <- raw_eRegulon_meta.df %>% dplyr::filter(class == "+/+")

# RSS BA data
RSS_mat.l <- readRDS(paste0(DATA_DIR, "RSS_BA_mats.rds"))
RSS_BA.m <- RSS_mat.l$RSS_BA.m
RSS_BA.z <- RSS_mat.l$RSS_BA.z

############################
### Process data & Setup ###
############################

# Get eRegulon & TFs
# Bascially  
N_TOP=50
top_eReg.v <- c()
for (i in 1:nrow(RSS_BA.z)) {
  top_eReg.v <- 
    colnames(RSS_BA.z)[order(RSS_BA.z[i,], decreasing = T)][1:N_TOP] %>%
    c(top_eReg.v, .) 
}
top_eReg.v <- unique(top_eReg.v)
if (length(top_eReg.v) > ncol(RSS_BA.z)) {top_eReg.v <- top_eReg.v[-length(top_eReg.v)]}
RSS_BA_top_Reg <- RSS_BA.z[, top_eReg.v]
colnames(RSS_BA_top_Reg) <- strsplit(top_eReg.v, '_') %>% 
  lapply(function(x) return(x[1])) %>%
  unlist()
## Set names
TOP_TF <- colnames(RSS_BA_top_Reg)
TOP_GENEeREG <- top_eReg.v
TOP_REGIONeREG <- BA_eRegulon_meta.df %>%
  filter(Gene_signature_name %in% top_eReg.v) %>%
  mutate(Gene_signature_name = factor(Gene_signature_name, levels = top_eReg.v)) %>%
  arrange(Gene_signature_name) %>%
  pull(Region_signature_name) %>%
  unique()

# Set TF expression matrix
seurat_TOP_TF <- subset(FL.SeuratObj, features = TOP_TF)
TF_EXPR.m <- GetAssayData(seurat_TOP_TF, layer = "data")

# Load AUCell scores
## AUC score gene
library(hdf5r)
h5file <- H5File$new("~/local_data/proj/Dev_Multiome/04.regulome/scp_ALL_pcw/lumi_outs/tmp_data/AllRegion_DAR_direct_gene_based_AUC.h5", mode = "r")
# h5file <- H5File$new(paste0(DATA_DIR, "AllRegion_DAR_direct_gene_based_AUC.h5"), mode = "r")
group <- h5file[["direct_gene_based_AUC"]]
row_names <- group[["axis0"]]$read()
col_names <- group[["axis1"]]$read()
auc_mat <- group[["block0_values"]]$read()
dimnames(auc_mat) <- list(row_names, col_names)
h5file2 <- H5File$new("~/local_data/proj/Dev_Multiome/04.regulome/scp_ALL_pcw/lumi_outs/tmp_data/AllRegion_DAR_extended_gene_based_AUC.h5", mode = "r")
# h5file <- H5File$new(paste0(DATA_DIR, "AllRegion_DAR_extended_gene_based_AUC.h5"), mode = "r")
group2 <- h5file2[["extended_gene_based_AUC"]]
row_names2 <- group2[["axis0"]]$read()
col_names2 <- group2[["axis1"]]$read()
auc_mat2 <- group2[["block0_values"]]$read()
dimnames(auc_mat2) <- list(row_names2, col_names2)
AUC_gene.m <- rbind(auc_mat, auc_mat2)[c(TOP_GENEeREG) , rownames(cell_metadata)]
### Set concensus names
AUC_GENE.m <- AUC_gene.m
rownames(AUC_GENE.m) <- TOP_TF
## AUC score region
h5file <- H5File$new("~/local_data/proj/Dev_Multiome/04.regulome/scp_ALL_pcw/lumi_outs/tmp_data/AllRegion_DAR_direct_region_based_AUC.h5", mode = "r")
# h5file <- H5File$new(paste0(DATA_DIR, "AllRegion_DAR_direct_region_based_AUC.h5"), mode = "r")
group <- h5file[["direct_region_based_AUC"]]
row_names <- group[["axis0"]]$read()
col_names <- group[["axis1"]]$read()
auc_mat <- group[["block0_values"]]$read()
dimnames(auc_mat) <- list(row_names, col_names)
h5file2 <- H5File$new("~/local_data/proj/Dev_Multiome/04.regulome/scp_ALL_pcw/lumi_outs/tmp_data/AllRegion_DAR_extended_region_based_AUC.h5", mode = "r")
# h5file <- H5File$new(paste0(DATA_DIR, "AllRegion_DAR_extended_region_based_AUC.h5"), mode = "r")
group2 <- h5file2[["extended_region_based_AUC"]]
row_names2 <- group2[["axis0"]]$read()
col_names2 <- group2[["axis1"]]$read()
auc_mat2 <- group2[["block0_values"]]$read()
dimnames(auc_mat2) <- list(row_names2, col_names2)
AUC_region.m <- rbind(auc_mat, auc_mat2)[c(TOP_REGIONeREG) , rownames(cell_metadata)]
### Set concensus names
AUC_REGION.m <- AUC_region.m
rownames(AUC_REGION.m) <- TOP_TF

# Save Data
saveRDS(list(TF_EXPR.m = TF_EXPR.m,
             AUC_gene.m = AUC_gene.m,
             AUC_GENE.m = AUC_GENE.m,
             AUC_region.m = AUC_region.m,
             AUC_REGION.m = AUC_REGION.m), 
        file = paste0(RES_DIR, "eRegulons_AUCell/TOPeReg_TF_AUCell_mat.rds"))

##################################################
### Training XGBoost for CT TFs identification ###
##################################################

# Function setting
source(paste0(WORK_DIR, "00.XGBoost_RMT_fun.R"))

# # Get RSS candidate eRegulons
# N_TOP = 80
# top_eRegulon.l <- list()
# for (i in CT_ORDER) {
#   top_eRegulon.l[[i]] <-
#     colnames(RSS_BA)[order(RSS_BA[i,], decreasing = T)][1:N_TOP] %>% 
#     strsplit(., '_') %>% 
#     lapply(function(x) return(x[1])) %>%
#     unlist() %>% 
#     intersect(., rownames(TF_EXPR.m))
# }

# Train XGBoost Model & Selection
train_TF_target.o <- train_OVR_models_targetCT(t(TF_EXPR.m),
                                               # unique_cts = "MK",
                                               cell_types = cell_metadata$anno_wnn_v51,
                                               model_type = "xgboost",
                                               nfold = 5,
                                               nrounds = 500,
                                               max_depth = 3,
                                               shap_ratio_cutoff = 2,
                                               shap_ratio_diff_cutoff = -Inf,
                                               seed = 123,
                                               do_zscore = TRUE,
                                               do_MinMax = FALSE,
                                               balance_classes = TRUE,
                                               n_iter = 30,
                                               parallel = TRUE,
                                               n_workers = 3, 
                                               nthread_xgb = 15)
# train_AUC_GENE_target.o <- train_OVR_models_targetCT(t(AUC_GENE.m),
#                                                      cell_types = cell_metadata$anno_wnn_v51,
#                                                      model_type = "xgboost",
#                                                      nfold = 5,
#                                                      nrounds = 500,
#                                                      max_depth = 3,
#                                                      shap_ratio_cutoff = 2,
#                                                      shap_ratio_diff_cutoff = -Inf,
#                                                      seed = 123,
#                                                      do_zscore = TRUE,
#                                                      do_MinMax = TRUE,
#                                                      balance_classes = TRUE,
#                                                      n_iter = 30,
#                                                      parallel = TRUE,
#                                                      n_workers = 3, 
#                                                      nthread_xgb = 15)
# train_AUC_REGION_target.o <- train_OVR_models_targetCT(t(AUC_REGION.m),
#                                                        cell_types = cell_metadata$anno_wnn_v51,
#                                                        model_type = "xgboost",
#                                                        nfold = 5,
#                                                        nrounds = 500,
#                                                        max_depth = 3,
#                                                        shap_ratio_cutoff = 2,
#                                                        shap_ratio_diff_cutoff = -Inf,
#                                                        seed = 123,
#                                                        do_zscore = TRUE,
#                                                        do_MinMax = TRUE,
#                                                        balance_classes = TRUE,
#                                                        n_iter = 30,
#                                                        parallel = TRUE,
#                                                        n_workers = 3, 
#                                                        nthread_xgb = 15)


# Save Res
# saveRDS(list(train_TF_target.o = train_TF_target.o,
#              train_AUC_GENE_target.o = train_AUC_GENE_target.o,
#              train_AUC_REGION_target.o = train_AUC_REGION_target.o),
#         file = paste0(RES_DIR, "eRegulons_CT_Filter/XGBoost_RMT_Res.rds"))
saveRDS(train_TF_target.o,
        file = paste0(RES_DIR, "eRegulons_CT_Filter/XGBoost_RMT_TFRes.rds"))

# Update selected feature ranks
XGBoost_RMT_TFRes.l <- readRDS(paste0(RES_DIR, "eRegulons_CT_Filter/XGBoost_RMT_TFRes.rds"))
for (i in 1:length(XGBoost_RMT_TFRes.l)) {
  tmp <- XGBoost_RMT_TFRes.l[[i]]
  shap.df <- tmp$shap
  ranked_features <- shap.df[order(-shap.df$SHAP_ratio), ]$Feature 
  tmp$ranked_features <- intersect(ranked_features, tmp$confident_features)
  XGBoost_RMT_TFRes.l[[i]] <- tmp
}
saveRDS(XGBoost_RMT_TFRes.l,
        file = paste0(RES_DIR, "eRegulons_CT_Filter/XGBoost_RMT_TFRes_wRankedFeature.rds"))
