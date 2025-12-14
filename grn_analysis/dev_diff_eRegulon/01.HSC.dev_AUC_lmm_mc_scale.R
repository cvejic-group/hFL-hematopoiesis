
# Initialization
setwd("/work/Local_Data/proj/Dev_Multiome/04.regulome_R/01.SCENICplus/05.SCENICplus_DevReg")
source("./00.LMM.R")
source("./00.Initialization.R")


library(tidyverse)
library(Seurat)
library(Matrix)
HSC.mc <- read.csv("/work/DevM_analysis/data/HSC_dev_diff_peaks/hsc_mc_bySample.csv")

dir.create(here::here("LMM_output_revisit"), showWarnings = FALSE)

# Load FL RNA data
FL.SeuratObj <- readRDS("~/local_data/proj/Dev_Multiome/data/FL_scrna_seurat_20250401.rds")
## Add WNN UMAP
FL_wnn.df <- read.csv("/work/DevM_analysis/01.annotation/10.integration_joint_clean/data/FL_wnn_umap.csv", 
                      row.names = 1)
colnames(FL.SeuratObj) <- gsub(pattern = "#", replacement = "_", 
                               x = colnames(FL.SeuratObj), fixed = TRUE)
FL_wnn.df <- FL_wnn.df[colnames(FL.SeuratObj) ,]
wnn_reduction <- CreateDimReducObject(embeddings = as.matrix(FL_wnn.df),
                                      key = "wnn_",
                                      assay = DefaultAssay(FL.SeuratObj))
FL.SeuratObj[["wnn"]] <- wnn_reduction
## Subset
FL_HSC.SeuratObj <- subset(FL.SeuratObj, cell = HSC.mc$barcode) |> 
  subset(subset = PCW != "Mixed")
rownames(HSC.mc) <- HSC.mc$barcode
FL_HSC.SeuratObj <- AddMetaData(FL_HSC.SeuratObj, HSC.mc)

# load mdata
mdata <- FL_HSC.SeuratObj@meta.data |>
  mutate(PCWsca = as.numeric(scale(as.integer(as.character(PCW))))) |>
  droplevels() |>
  as.data.frame()
str(mdata)

# Load AUCell selection
RSS_BA_mats <- readRDS("~/local_data/proj/Dev_Multiome/04.regulome_R/04.SCENICplus_RawSummary/results/RSS_BA_mats.rds")
top_eReg.v <- colnames(RSS_BA_mats$RSS_BA.m)
top_tf <- strsplit(top_eReg.v, split = "_") |> 
  lapply(function(x) return(x[1])) |> 
  unlist()

# Load AUC data
library(hdf5r)
h5file <- H5File$new("~/local_data/proj/Dev_Multiome/04.regulome/scp_ALL_pcw/lumi_outs/tmp_data/AllRegion_DAR_direct_gene_based_AUC.h5", mode = "r")
group <- h5file[["direct_gene_based_AUC"]]
row_names <- group[["axis0"]]$read()
col_names <- group[["axis1"]]$read()
auc_mat <- group[["block0_values"]]$read()
dimnames(auc_mat) <- list(row_names, col_names)
h5file2 <- H5File$new("~/local_data/proj/Dev_Multiome/04.regulome/scp_ALL_pcw/lumi_outs/tmp_data/AllRegion_DAR_extended_gene_based_AUC.h5", mode = "r")
group2 <- h5file2[["extended_gene_based_AUC"]]
row_names2 <- group2[["axis0"]]$read()
col_names2 <- group2[["axis1"]]$read()
auc_mat2 <- group2[["block0_values"]]$read()
dimnames(auc_mat2) <- list(row_names2, col_names2)
AUC_gene.m <- rbind(auc_mat, auc_mat2)[c(top_eReg.v) ,]
## Change cell name
colnames(AUC_gene.m) <- gsub(pattern = "#", replacement = "_", 
                               x = colnames(AUC_gene.m), fixed = TRUE)
AUC_gene_HSC.m <- AUC_gene.m[, rownames(mdata)]
rownames(AUC_gene_HSC.m) <- top_tf

colnames(AUC_region.m) <- gsub(pattern = "#", replacement = "_", 
                             x = colnames(AUC_region.m), fixed = TRUE)
AUC_region_HSC.m <- AUC_region.m[, rownames(mdata)]
rownames(AUC_region_HSC.m) <- top_tf

# peak mat
library(Matrix)
agg_counts <- t(sapply(split(seq_along(mdata$mc_group), mdata$mc_group), 
                       function(idx) {
  Matrix::rowSums(AUC_gene_HSC.m[, idx, drop = FALSE])
}))
agg_counts <- t(agg_counts)

Y_scale <- t(scale(t(AUC_gene_HSC.m)))
stopifnot(all(colnames(Y_scale) == mdata$barcode))

# linear mixed model
mdata$log_nCount_RNA <- log1p(mdata$nCount_RNA)
mdata$log_nFeature_RNA <- log1p(mdata$nFeature_RNA)
res <- LFLMM(Y_scale, mdata[,c("log_nCount_RNA", "log_nFeature_RNA", "percent.mt", "libraryID", "donorID",
                               "Sex", "Batch", "PCWsca")], ITRMAX=300)
saveRDS(res, file = here::here("LMM_output_revisit", "HSC_AUC_region_lmm_res.rds"))

# de with scale.data
DOCNAME <- ""
if (TRUE) {
  # de
  de <- getBF(Y_scale, res, "PCWsca", DE1 = NA)
  names(de)
  # de df
  df_de <- data.frame(gene = rownames(de$beta), beta = de$beta[,1], ltsr = de$ltsr[,1]) |>
    arrange(-beta)
  saveRDS(de, file = here::here("LMM_output_revisit", DOCNAME, "HSC_dev_AUC_gene_lmm.rds"))
  saveRDS(df_de, file = here::here("LMM_output_revisit", DOCNAME, "HSC_dev_AUC_gene_lmm_df.rds"))
}


# df_de2 <- df_de %>% dplyr::filter(ltsr >= 0.8)
# 
# HSC_lmm_de_time_ltsr0_9 %>% dplyr::filter(gene %in% rownames(df_de2))
# plot(HSC_lmm_de_time_df[rownames(df_de2),"beta"], df_de2$beta)
# abline(0,1)
# abline(0,1e100)
# abline(0,0)
# cor.test(HSC_lmm_de_time_df[rownames(df_de2),"beta"], df_de2$beta)

df_de <- readRDS("./LMM_output/HSC_dev_AUC_region_lmm_df.rds")
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
  # xlim(-.5,.5) +
  theme_classic() +
  guides(color = "none",#color = guide_legend(override.aes = list(size = 5)),
         size = "none") +
  theme(
    axis.title = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", hjust = 0, size = 11),
    plot.title = element_text(face = "bold", hjust = .5, size = 15)
  )
dev.off()

HSC_lmm_de_time_df <- readRDS("/work/DevM_analysis/data/HSC_dev_diff_genes/HSC_lmm_de_time_df.rds")
HSC_lmm_de_time_df2 <- HSC_lmm_de_time_df[df_de$gene ,]
df_de_gene <- HSC_lmm_de_time_df2#df_de_gene[df_de_region$gene ,]


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
    title = "Gene-based eRegulon Dynamics Coordinates with TF Dynamics",
    x = "Beta_eRegulon_GeneBased",
    y = "Beta_Gene",
    size = "Average LTSR"
  ) +
  # xlim(-0.55, .65) +
  # ylim(-0.2, .25) +
  theme_classic() +
  guides(color = "none",#color = guide_legend(override.aes = list(size = 5)),
         )+#size = "none") +
  theme(
    axis.title = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", hjust = 0, size = 11),
    plot.title = element_text(face = "bold", hjust = .5, size = 15)
  )
dev.off()

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
  # xlim(-0.55, .65) +
  # ylim(-0.2, .25) +
  theme_classic() +
  guides(color = "none",#color = guide_legend(override.aes = list(size = 5)),
  )+#size = "none") +
  theme(
    axis.title = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", hjust = 0, size = 11),
    plot.title = element_text(face = "bold", hjust = .5, size = 15)
  )
dev.off()


# all_df <- all_df %>% 
#   dplyr::filter(ltsr >= 0.5 | ltsr_gene >= 0.5)
df_long <- all_df %>% 
  pivot_longer(cols = c(beta, beta_gene), names_to = "type", values_to = "value") %>%
  mutate(ltsr_val = ifelse(type == "beta", ltsr, ltsr_gene),
         type = recode(type, beta = "Region-level", beta_gene = "Gene-level"))
all_df <- all_df[] %>%
  mutate(sim_score = abs(beta - beta_gene))

gene_order <- all_df %>%
  arrange(sim_score) %>%
  pull(gene)

df_long$gene <- factor(df_long$gene, levels = gene_order)
Cairo::CairoPDF("eRegulon_Stat_Comparison.pdf", width = 28, height = 2.5, family = "Arial")
ggplot(df_long, aes(x = gene, y = type)) +
  geom_tile(aes(fill = value), color = "white") +
  geom_point(aes(size = ltsr_val), shape = 21, fill = "black") +
  scale_fill_gradientn(
    colours = c("#3a86ff", "#48cae4", "white", "#e9c46a", "#d62828"),
    values = scales::rescale(c(0, .25, .5, .75, 1)),
    name = "Effect size"
  ) +
  scale_size(range = c(1, 4), name = "LTSR") +
  scale_x_discrete(limits = rev(gene_order)) +
  theme_minimal(base_size = 12) +
  labs(x = NULL, y = NULL, title = "Comparison of Effect size and LTSR") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "top",
        axis.text = element_text(face = "bold"),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = .5, colour = "black", size = 14))
  # coord_flip()
dev.off()




# # Try CCAT
# library(clusterProfiler)
# library(org.Hs.eg.db)
# gene_entrez <- bitr(unique(rownames(FL_HSC.SeuratObj)), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
# data.mat <- GetAssayData(FL_HSC.SeuratObj, layer = "counts")
# data.mat <- data.mat[gene_entrez$SYMBOL ,]
# rownames(data.mat) <- gene_entrez$ENTREZID
# ccat.o <- CCAT(data.m = data.mat, ppiA.m = net13Jun12.m, log_trans = T, parallelMode = T, mcores = 5, subsamplesize = 1000)
# head(ccat.o)
# FL_HSC.SeuratObj$CCAT <- ccat.o
# ggplot(FL_HSC.SeuratObj@meta.data, aes(Batch, CCAT)) + geom_boxplot()# + scale_x_discrete(limits = as.character(5:18))


