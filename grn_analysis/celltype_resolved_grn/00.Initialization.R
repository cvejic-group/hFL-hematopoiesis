#################################################
### SCENIC+ results: Cell Type Specific eGRNs ###
#################################################

###-------------------Note-------------------###
###-------      Run with R-4.4.2      -------###
###------------------------------------------###

##################
# Initialization #
##################

# Set PATH
setwd("~/local_data/proj/Dev_Multiome/04.regulome_R/01.SCENICplus/04.SCENICplus_CTeGRN/")
WORK_DIR = "~/local_data/proj/Dev_Multiome/04.regulome_R/01.SCENICplus/04.SCENICplus_CTeGRN//"
DATA_DIR = paste0(WORK_DIR, "data/")
RES_DIR = paste0(WORK_DIR, "results/")
FIG_DIR = paste0(WORK_DIR, "plots/")
SCP_RES="~/local_data/proj/Dev_Multiome/04.regulome/scp_ALL_pcw/lumi_outs/tmp_data/"

# Load essential packages
library(tidyverse)
library(tidyr)
library(tibble)
library(data.table)
library(rhdf5)
library(Seurat)
library(SeuratDisk)
library(SingleCellExperiment)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(patchwork)

# Cell type order
CT_ORDER=c("HSC", "GP", "Granulocyte",
           "MEMP-t", "MEMP", "MEP", "MEMP-Mast-Ery", "MEMP-Ery", "Early-Ery", "Late-Ery",
           "MEMP-MK", "MK", "MastP-t", "MastP", "Mast",
           "MDP", "Monocyte", "Kupffer", "cDC1", "cDC2", "pDC", "ASDC",
           "LMPP", "LP", "Cycling-LP", "PreProB", "ProB-1", "ProB-2", "Large-PreB", "Small-PreB", "IM-B",
           "NK", "ILCP", "T",
           "Hepatocyte", "Endothelia")
## Set COLs
COL <- c("#E41A1C", "#E0FFFF", "#B3CDE3", "#E6AB02", "#FF7F00",
         "#CD661D", "#FDCDAC", "#E9967A", "#CD5555", "#8B0000",
         "#663C1F", "#40E0D0", "#1E90FF", "#1F78B4", "#253494",
         "#E6F5C9", "#005A32", "#00EE00", "#ADFF2F", "#B3DE69",
         "#4DAF4A", "#CDC673", "#FFF2AE", "#FFD92F", "#FFFF33",
         "#FFF0F5", "#FFB5C5", "#E78AC3", "#CD1076", "#FF3E96",
         "#FF00FF", "#A020F0", "#49006A", "#984EA3", "#666666",
         "#000000")
names(COL) <- CT_ORDER
COL4 <- c("#c1121f", "#ffb703", "#669bbc", "#606c38")

