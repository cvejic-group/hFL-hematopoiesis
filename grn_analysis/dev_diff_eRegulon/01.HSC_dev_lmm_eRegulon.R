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
setwd("~/work/")
source("./00.LMM.R")
source("./00.Initialization.R")
library(tidyverse)
library(Seurat)
library(Matrix)

dir.create(here::here("LMM_output"), showWarnings = FALSE)

# Load FL RNA data
FL.SeuratObj <- readRDS("~/work/FL_scrna_seurat.rds")
## Add WNN UMAP
FL_wnn.df <- read.csv("FL_wnn_umap.csv", 
                      row.names = 1)
colnames(FL.SeuratObj) <- gsub(pattern = "#", replacement = "_", 
                               x = colnames(FL.SeuratObj), fixed = TRUE)
FL_wnn.df <- FL_wnn.df[colnames(FL.SeuratObj) ,]
wnn_reduction <- CreateDimReducObject(embeddings = as.matrix(FL_wnn.df),
                                      key = "wnn_",
                                      assay = DefaultAssay(FL.SeuratObj))
FL.SeuratObj[["wnn"]] <- wnn_reduction
## Subset
FL_HSC.SeuratObj <- subset(FL.SeuratObj, cell = celltype == "HSC") |> 
  subset(subset = PCW != "Mixed")
mdata <- FL_HSC.SeuratObj@meta.data |>
  mutate(PCWsca = as.numeric(scale(as.integer(as.character(PCW))))) |>
  droplevels() |>
  as.data.frame()

# Load AUC data
# Load data
TOPeReg_TF_AUCell_mat <- readRDS(paste0(RES_DIR, "TOPeReg_TF_AUCell_mat.rds"))
TF_EXPR.m <- TOPeReg_TF_AUCell_mat$TF_EXPR.m
AUC_GENE.m <- TOPeReg_TF_AUCell_mat$AUC_GENE.m
AUC_REGION.m <- TOPeReg_TF_AUCell_mat$AUC_REGION.m
## Change cell name
colnames(AUC_GENE.m) <- gsub(pattern = "#", replacement = "_", 
                             x = colnames(AUC_GENE.m), fixed = TRUE)
AUC_gene_HSC.m <- AUC_GENE.m[, rownames(mdata)]
colnames(AUC_REGION.m) <- gsub(pattern = "#", replacement = "_", 
                               x = colnames(AUC_REGION.m), fixed = TRUE)
AUC_region_HSC.m <- AUC_REGION.m[, rownames(mdata)]

# Scale AUCell scores
Y_scale <- t(scale(t(AUC_gene_HSC.m)))
stopifnot(all(colnames(Y_scale) == mdata$barcode))

# linear mixed model
mdata$log_nCount_RNA <- log1p(mdata$nCount_RNA)
mdata$log_nFeature_RNA <- log1p(mdata$nFeature_RNA)
res <- LFLMM(Y_scale, mdata[,c("log_nCount_RNA", "log_nFeature_RNA", "percent.mt", "libraryID", "donorID",
                               "Sex", "Batch", "PCWsca")], ITRMAX=300)
saveRDS(res, file = here::here("LMM_output", "HSC_AUC_region_lmm_res.rds"))

# de with scale.data
DOCNAME <- ""
if (TRUE) {
  # de
  de <- getBF(Y_scale, res, "PCWsca", DE1 = NA)
  names(de)
  # de df
  df_de <- data.frame(gene = rownames(de$beta), beta = de$beta[,1], ltsr = de$ltsr[,1]) |>
    arrange(-beta)
  saveRDS(de, file = here::here("LMM_output", DOCNAME, "HSC_dev_AUC_gene_lmm.rds"))
  saveRDS(df_de, file = here::here("LMM_output", DOCNAME, "HSC_dev_AUC_gene_lmm_df.rds"))
}

# Generate plots
## Set threshold
df_de$significant <- df_de$ltsr > 0.75
N_top <- length(which(df_de$significant))
top_select <- df_de[order(df_de$ltsr, decreasing = T), ][1:N_top, ]
Cairo::CairoPDF("SuppFig16a.pdf", width = 12, height = 12, family = "Arial")
ggplot(df_de, aes(x = beta, y = ltsr, color = significant, size = ltsr)) +
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

# Coordination with gene's LMM
HSC_lmm_de_time_df <- readRDS("~/work/HSC_lmm_de_time_df.rds")
HSC_lmm_de_time_df2 <- HSC_lmm_de_time_df[df_de$gene ,]
df_de_gene <- HSC_lmm_de_time_df2
TF.idx <- intersect(rownames(df_de), rownames(df_de_gene))

all_df <- df_de[TF.idx ,]
all_df$beta_gene <- df_de_gene[TF.idx ,]$beta
all_df$ltsr_avg <- (all_df$ltsr + df_de_gene[TF.idx ,]$ltsr)/2
all_df$significant <- all_df$ltsr > 0.75
N_top <- length(which(all_df$significant))
top_select <- all_df[order(all_df$ltsr, decreasing = T), ][1:N_top, ]
Cairo::CairoPDF("SuppFig16b.pdf", width = 12, height = 12, family = "Arial")
ggplot(all_df, aes(x = beta, y = beta_gene, color = significant, size = ltsr_avg)) +
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

# Coordination with region_eRegulon's LMM
df_de_gene <- readRDS("./LMM_output/HSC_dev_AUC_lmm_df.rds")
all_df <- df_de
all_df$beta_gene <- df_de_gene[all_df$gene,]$beta
all_df$ltsr_gene <- df_de_gene[all_df$gene,]$ltsr

all_df$ltsr_avg <- (all_df$ltsr + all_df$ltsr_gene)/2
all_df$significant <- (all_df$ltsr > 0.75 | all_df$ltsr_gene > 0.75)
N_top <- length(which(all_df$significant))
top_select <- all_df[which(all_df$significant), ]
Cairo::CairoPDF("SuppFig16c.pdf", width = 12, height = 12, family = "Arial")
ggplot(all_df, aes(x = beta, y = beta_gene, color = significant, size = ltsr_avg)) +
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
