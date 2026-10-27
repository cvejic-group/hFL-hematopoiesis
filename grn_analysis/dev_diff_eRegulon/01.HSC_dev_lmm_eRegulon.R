############################################
### eRegulon activity across development ###
############################################

###-------------------Note-------------------###
###-------      Run with R-4.4.2      -------###
###------------------------------------------###

##################
# Initialization #
##################

# Initialization
setwd("/work/")
source("./00.LMM.R")
source("./00.Initialization.R")
library(tidyverse)
library(Seurat)
library(Matrix)

dir.create(here::here("LMM_output"), showWarnings = FALSE)
DOCNAME <- ""

#############################
# Part 1: eREgulon Analysis #
#############################

# Load FL RNA data
FL.SeuratObj <- readRDS("~/work/FL_scrna_seurat.rds")
## Subset
FL_HSC.SeuratObj <- subset(FL.SeuratObj, cell = celltype == "HSC") |> 
  subset(subset = PCW != "Mixed")
mdata <- FL_HSC.SeuratObj@meta.data |>
  mutate(PCWsca = as.numeric(scale(as.integer(as.character(PCW))))) |>
  droplevels() |>
  as.data.frame()

# Load AUC data
# Load data
TOPeReg_TF_AUCell_mat <- readRDS("TOPeReg_TF_AUCell_mat.rds")
TF_EXPR.m <- TOPeReg_TF_AUCell_mat$TF_EXPR.m
AUC_GENE.m <- TOPeReg_TF_AUCell_mat$AUC_GENE.m
AUC_REGION.m <- TOPeReg_TF_AUCell_mat$AUC_REGION.m
## Change cell name if needed
# colnames(AUC_GENE.m) <- gsub(pattern = "#", replacement = "_", 
#                              x = colnames(AUC_GENE.m), fixed = TRUE)
AUC_gene_HSC.m <- AUC_GENE.m[, rownames(mdata)]
# colnames(AUC_REGION.m) <- gsub(pattern = "#", replacement = "_", 
#                                x = colnames(AUC_REGION.m), fixed = TRUE)
AUC_region_HSC.m <- AUC_REGION.m[, rownames(mdata)]

# Scale AUCell scores
Y_eReg_gene_scale <- t(scale(t(AUC_gene_HSC.m)))
stopifnot(all(colnames(Y_eReg_gene_scale) == mdata$barcode))
Y_eReg_region_scale <- t(scale(t(AUC_region_HSC.m)))
stopifnot(all(colnames(Y_eReg_region_scale) == mdata$barcode))

# Gene-based eRegulon: linear mixed model
res_eReg_gene <- LFLMM(Y_eReg_gene_scale, mdata[,c("percent.mt", "libraryID", "donorID",
                                                   "Sex", "Batch", "PCWsca")], ITRMAX=300)
saveRDS(res_eReg_gene, file = here::here("LMM_output", "HSC_AUC_gene_lmm_res.rds"))

# Gene-based eRegulon: de with scale.data
if (TRUE) {
  # de
  de_eReg_gene <- getBF(Y_eReg_gene_scale, res_eReg_gene, "PCWsca", DE1 = NA)
  names(de_eReg_gene)
  # de df
  df_eReg_gene <- data.frame(gene = rownames(de_eReg_gene$beta), beta = de_eReg_gene$beta[,1], ltsr = de_eReg_gene$ltsr[,1]) |>
    arrange(-beta)
  saveRDS(de_eReg_gene, file = here::here("LMM_output", DOCNAME, "HSC_dev_AUC_gene_lmm.rds"))
  saveRDS(df_eReg_gene, file = here::here("LMM_output", DOCNAME, "HSC_dev_AUC_gene_lmm_df.rds"))
}

# Region-based eRegulon: linear mixed model
res_eReg_region <- LFLMM(Y_eReg_region_scale, mdata[,c("percent.mt", "libraryID", "donorID",
                                                       "Sex", "Batch", "PCWsca")], ITRMAX=300)
saveRDS(res_eReg_region, file = here::here("LMM_output", "HSC_AUC_region_lmm_res.rds"))

# Region-based eRegulon: de with scale.data
if (TRUE) {
  # de
  de_eReg_region <- getBF(Y_eReg_region_scale, res_eReg_region, "PCWsca", DE1 = NA)
  names(de_eReg_region)
  # de df
  df_eReg_region <- data.frame(gene = rownames(de_eReg_region$beta), beta = de_eReg_region$beta[,1], ltsr = de_eReg_region$ltsr[,1]) |>
    arrange(-beta)
  saveRDS(de_eReg_region, file = here::here("LMM_output", DOCNAME, "HSC_dev_AUC_region_lmm.rds"))
  saveRDS(df_eReg_region, file = here::here("LMM_output", DOCNAME, "HSC_dev_AUC_region_lmm_df.rds"))
}

####################
# Part 2: Plotting #
####################

# Load LMM results
df_eReg_gene <- readRDS(here::here("LMM_output", DOCNAME, "HSC_dev_AUC_gene_lmm_df.rds"))
df_eReg_region <- readRDS(here::here("LMM_output", DOCNAME, "HSC_dev_AUC_region_lmm_df.rds"))
df_TF_expr <- readRDS("~/work/HSC_lmm_de_time_df.rds")

# Region-based eRegulon dynamics across development time
## Set threshold
df_eReg_region$significant <- df_eReg_region$ltsr > 0.75
N_top <- length(which(df_eReg_region$significant))
top_select <- df_eReg_region[order(df_eReg_region$ltsr, decreasing = T), ][1:N_top, ]
Cairo::CairoPDF("SuppFig_a.pdf", width = 12, height = 12, family = "Arial")
ggplot(df_eReg_region, aes(x = beta, y = ltsr, color = significant, size = ltsr)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("grey70", "#253494"),
                     name = "LTSR significant") +
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "grey30") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
  geom_text_repel(
    data = top_select,
    aes(label = gene),
    size = 3,
    max.overlaps = Inf,
    min.segment.length = 0.1,
    box.padding = .7
  ) +
  labs(
    title = "Region-based eRegulon Dynamics across Development Time",
    x = "Beta",
    y = "LTSR"
  ) +
  ylim(0, 1.05) +
  theme_classic() +
  guides(color = "none",
         size = "none") +
  theme(
    axis.title = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", hjust = 0, size = 11),
    plot.title = element_text(face = "bold", hjust = .5, size = 15)
  )
dev.off()

# Coordination of gene-based eRegulon with TF expression LMM
df_TF_expr <- df_TF_expr[df_eReg_gene$gene ,]
TF.idx <- intersect(rownames(df_eReg_gene), rownames(df_TF_expr))

df_eReg_gene_vs_TF <- df_eReg_gene[TF.idx ,]
df_eReg_gene_vs_TF$beta_TF <- df_TF_expr[TF.idx ,]$beta
df_eReg_gene_vs_TF$ltsr_avg <- (df_eReg_gene_vs_TF$ltsr + df_TF_expr[TF.idx ,]$ltsr)/2
df_eReg_gene_vs_TF$significant <- df_eReg_gene_vs_TF$ltsr > 0.75
N_top <- length(which(df_eReg_gene_vs_TF$significant))
top_select <- df_eReg_gene_vs_TF[order(df_eReg_gene_vs_TF$ltsr, decreasing = T), ][1:N_top, ]
Cairo::CairoPDF("Fig_e.pdf", width = 12, height = 12, family = "Arial")
ggplot(df_eReg_gene_vs_TF, aes(x = beta, y = beta_TF, color = significant, size = ltsr_avg)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("grey70", "#253494"),
                     name = "LTSR significant") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_text_repel(
    data = top_select,
    aes(label = gene),
    size = 3,
    max.overlaps = Inf,
    min.segment.length = 0.1,
    box.padding = .7
  ) +
  labs(
    title = "Gene-based eRegulon Dynamics Coordinates with TF Dynamics",
    x = "Beta_eRegulon_GeneBased",
    y = "Beta_Gene",
    size = "Average LTSR"
  ) +
  theme_classic() +
  guides(color = "none") +
  theme(
    axis.title = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", hjust = 0, size = 11),
    plot.title = element_text(face = "bold", hjust = .5, size = 15)
  )
dev.off()

# Coordination of region-based eRegulon with gene-based eRegulon LMM
df_eReg_region_vs_gene <- df_eReg_region
df_eReg_region_vs_gene$beta_gene <- df_eReg_gene[df_eReg_region_vs_gene$gene,]$beta
df_eReg_region_vs_gene$ltsr_gene <- df_eReg_gene[df_eReg_region_vs_gene$gene,]$ltsr

df_eReg_region_vs_gene$ltsr_avg <- (df_eReg_region_vs_gene$ltsr + df_eReg_region_vs_gene$ltsr_gene)/2
df_eReg_region_vs_gene$significant <- (df_eReg_region_vs_gene$ltsr > 0.75 | df_eReg_region_vs_gene$ltsr_gene > 0.75)
N_top <- length(which(df_eReg_region_vs_gene$significant))
top_select <- df_eReg_region_vs_gene[which(df_eReg_region_vs_gene$significant), ]
Cairo::CairoPDF("SuppFig_c.pdf", width = 12, height = 12, family = "Arial")
ggplot(df_eReg_region_vs_gene, aes(x = beta, y = beta_gene, color = significant, size = ltsr_avg)) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = c("grey70", "#253494"),
                     name = "LTSR significant") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  # geom_abline(slope = 1) +
  geom_text_repel(
    data = top_select,
    aes(label = gene),
    size = 3,
    max.overlaps = Inf,
    min.segment.length = 0.1,
    box.padding = .7
  ) +
  labs(
    title = "Coordination of Two eRegulon Dynamics",
    x = "Beta_eRegulon_RegionBased",
    y = "Beta_eRegulon_GeneBased",
    size = "Average LTSR"
  ) +
  theme_classic() +
  guides(color = "none",
  ) +
  theme(
    axis.title = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", hjust = 0, size = 11),
    plot.title = element_text(face = "bold", hjust = .5, size = 15)
  )
dev.off()