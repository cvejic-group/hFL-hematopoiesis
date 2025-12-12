########################################################
###        Exteral Data - Calvanese 2022 Nature      ###
###        Label Transfer for our PCW5 HSC data      ###
########################################################

# Set working directory
setwd("~/local_data/proj/Dev_Multiome/05.External_Label_Transfer/")

# Load essential packages
library(scater)
library(Seurat)
library(Signac)
library(ggplot2)
library(cowplot)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

# Set COL
CT_ORDER=c("HSC", "GP", "Granulocyte",
           "MEMP-t", "MEMP", "MEP", "MEMP-Mast-Ery", "MEMP-Ery", "Early-Ery", "Late-Ery",
           "MEMP-MK", "MK", "MastP-t", "MastP", "Mast",
           "MDP", "Monocyte", "Kupffer", "cDC1", "cDC2", "pDC", "ASDC",
           "LMP", "LP", "Cycling-LP", "PreProB", "ProB-1", "ProB-2", "Large-PreB", "Small-PreB", "IM-B",
           "NK", "ILCP", "T",
           "Hepatocyte", "Endothelia")
COL <- c("#E41A1C",
         "#E0FFFF",
         "#B3CDE3",
         "#E6AB02",
         "#FF7F00",
         "#CD661D",
         "#FDCDAC",
         "#E9967A",
         "#CD5555",
         "#8B0000",
         "#663C1F",
         "#40E0D0",
         "#1E90FF",
         "#1F78B4",
         "#253494",
         "#E6F5C9",
         "#005A32",
         "#00EE00",
         "#ADFF2F",
         "#B3DE69",
         "#4DAF4A",
         "#CDC673",
         "#FFF2AE",
         "#FFD92F",
         "#FFFF33",
         "#FFF0F5",
         "#FFB5C5",
         "#E78AC3",
         "#CD1076",
         "#FF3E96",
         "#FF00FF",
         "#A020F0",
         "#49006A",
         "#984EA3",
         "#666666",
         "#000000")
names(COL) <- CT_ORDER
COL8 <- c("#E41A1C","#FF7F00","#CD1076", "#F781BF", "#984EA3", "#377EB8","#4DAF4A", "#A65628")
names(COL8) <- c("HSC",
                 "Lympho",
                 "Granulo",
                 "Gr/Prolif",
                 "Mo/Mϕ/Prolif",
                 "Mo/Mϕ",
                 "Mo/Mϕ/Endo",
                 "Mo/Mϕ/Stroma")

##############################
### Subset Process FL PCW5 ###
##############################

# Load FL data
FL.SeuratObj <- readRDS("~/local_data/proj/Dev_Multiome/data/FL_scrna_seurat_20250401.rds")
FL.SeuratObj$covar_group <- paste(FL.SeuratObj$donorID, FL.SeuratObj$libraryID, sep = "_")

# Subset to HSC
FL_PCW5_HSC.SeuratObj <- subset(FL.SeuratObj, subset = anno_wnn_v51 == "HSC")
## Normalization
FL_PCW5_HSC.SeuratObj <- NormalizeData(FL_PCW5_HSC.SeuratObj)

# Save HSC Object
saveRDS(FL_PCW5_HSC.SeuratObj, "./FL_PCW5_HSC_SeuratObj.rds")

############################################
### Process Calvanese et al. AGM CS14&15 ###
############################################

# Load Calvanese_EDFig1f_AGMHema data as reference
## EDFig-1f re-process
Calvanese_EDFig1f_AGMHema.SeuratObj <- readRDS("~/Local_Data/proj/Dev_Multiome/data/Calvanese_EDFig1f_SeuratObj.rds")
### Run UMAP with settings specified by original paper
Calvanese_EDFig1f_AGMHema.SeuratObj <- RunUMAP(Calvanese_EDFig1f_AGMHema.SeuratObj, 
                                               dims = 1:23, 
                                               reduction.name = "umap_rerun", 
                                               return.model = T)
## Re-name the seurat clusters based on EDFig1f
Calvanese_EDFig1f_AGMHema.SeuratObj@meta.data <- Calvanese_EDFig1f_AGMHema.SeuratObj@meta.data %>%
  mutate(CT = recode(seurat_clusters,
                     "1" = "HSC",
                     "7" = "Lympho",
                     "3" = "Granulo",
                     "5" = "Gr/Prolif",
                     "0" = "Mo/Mϕ/Prolif",
                     "2" = "Mo/Mϕ",
                     "4" = "Mo/Mϕ",
                     "6" = "Mo/Mϕ/Endo",
                     "8" = "Mo/Mϕ/Stroma"))
## Plot AGM Ori UMAP
cellNum.idx <- table(Calvanese_EDFig1f_AGMHema.SeuratObj$CT)
DimPlot(
  Calvanese_EDFig1f_AGMHema.SeuratObj,
  group.by = "CT",
  reduction = "umap_rerun",
  pt.size = 1, 
  alpha = 1,
  label = T,
  repel = T, 
  label.size = 5
) + 
  scale_color_manual(values = COL8, labels = paste0(names(cellNum.idx), " (", cellNum.idx, ")")) +
  labs(title = "Original UMAP of Calvanese et al.\n EDFig1f AGM Hematopoietic") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = .5, colour = "black", size = 17))

# Save data
saveRDS(Calvanese_EDFig1f_AGMHema.SeuratObj, "./Calvanese_EDFig1f_Reprocess_SeuratObj.rds")

####################################
### Construct symphony Reference ###
####################################

library(symphony)

# Build reference
set.seed(1)
AGM_CS14.15_sym.reference = symphony::buildReference(
  GetAssayData(Calvanese_EDFig1f_AGMHema.SeuratObj, layer = "counts"),
  Calvanese_EDFig1f_AGMHema.SeuratObj@meta.data,
  vars = NULL,                # variables to integrate over
  K = 100,                    # number of Harmony clusters
  verbose = TRUE,             # verbose output
  do_umap = TRUE,             # can set to FALSE if want to run umap separately later
  do_normalize = TRUE,        # set to TRUE if input counts are not normalized yet
  vargenes_method = 'vst',    # method for variable gene selection ('vst' or 'mvp')
  vargenes_groups = NULL,     # metadata column specifying groups for variable gene selection 
  topn = 2000,                # number of variable genes to choose per group
  d = 23,                     # number of PCs
  save_uwot_path = './sym_Ref_AGM_uwot_model_DR23'
)
AGM_CS14.15_sym.reference$normalization_method = 'log(CP10k+1)' # optionally save normalization method in custom slot
## Save reference
saveRDS(AGM_CS14.15_sym.reference, './sym_Ref_AGM_DR23.rds')

####################################################
## Plot Reference data with original annotations ###
####################################################

# Prepare plot data
umap_labels_SymRef.df <- cbind(Calvanese_EDFig1f_AGMHema.SeuratObj@meta.data, 
                               AGM_CS14.15_sym.reference$umap$embedding)
cellNum.idx <- table(umap_labels_SymRef.df$CT)
cellNum.idx <- cellNum.idx[names(COL8)]
# Set arrow
arrow.l <- list(x = -7, y = -8, x_len = 2, y_len = 1.5)
## Plot
AGM_p3.1 <- ggplot() +
  geom_point(data = umap_labels_SymRef.df, 
             aes(UMAP1, UMAP2, color = CT), alpha = 1, size = .1) +
  scale_color_manual(values = COL8[names(cellNum.idx)], 
                     breaks = names(cellNum.idx), 
                     name = "Calvanese et al. Annotation\n(cell number)",
                     labels = paste0(names(cellNum.idx), " (", cellNum.idx, ")")) +
  labs(title = "Symphony Constructed Reference\n(Calvanese et al. AGM CS14 & CS15)") +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  annotate("segment",
           x = arrow.l$x, xend = arrow.l$x + c(arrow.l$x_len, 0),
           y = arrow.l$y, yend = arrow.l$y + c(0, arrow.l$y_len),
           arrow = arrow(type = "open", length = unit(.1, 'inches')),
           linewidth = .5) +
  annotate("text",
           x = arrow.l$x + .9, y = arrow.l$y + .4, label = "UMAP", 
           size = 3.5,
           color = "black",
           fontface = "bold") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", hjust = .5, colour = "black", size = 13),
        legend.title = element_text(face = "bold", hjust = 0.5, size = 11),
        text = element_text(family = "Noto Sans"))
## Save plot
Cairo::CairoPDF("AGM_p3.1.pdf", width = 8, height = 6, family = "Arial")
AGM_p3.1
dev.off()

###############################
### symphony Label Transfer ###
###############################

# Map query
HSC2AGM_CS14.15.query = mapQuery(GetAssayData(FL_PCW5_HSC.SeuratObj, layer = "counts"),       # query gene expression (genes x cells)
                                 FL_PCW5_HSC.SeuratObj@meta.data,                             # query metadata (cells x attributes)
                                 AGM_CS14.15_sym.reference,                                   # Symphony reference object
                                 vars = 'covar_group',                                        # Correction for technical factors
                                 do_normalize = TRUE,                                         # perform log(CP10k+1) normalization on query
                                 do_umap = TRUE)                                              # project query cells into reference UMAP

# Predicted Labels
HSC2AGM_CS14.15.query <- knnPredict(HSC2AGM_CS14.15.query, AGM_CS14.15_sym.reference, 
                                    AGM_CS14.15_sym.reference$meta_data$CT, k = 7)
table(HSC2AGM_CS14.15.query$meta_data$cell_type_pred_knn)
## Update SeuratObject
FL_PCW5_HSC.SeuratObj$AGM_predicted_label <- as.character(HSC2AGM_CS14.15.query$meta_data$cell_type_pred_knn)
FL_PCW5_HSC.SeuratObj$AGM_predicted_scores <- as.numeric(HSC2AGM_CS14.15.query$meta_data$cell_type_pred_knn_prob)

#######################################
## Visualization of data projection ###
#######################################

# Prepare data
bk.umap.AGM <- cbind(Calvanese_EDFig1f_AGMHema.SeuratObj@meta.data, 
                     AGM_CS14.15_sym.reference$umap$embedding)
bk.umap.AGM$Label <- "Background"
FL2AGM.umap.sym <- as.data.frame(HSC2AGM_CS14.15.query$umap)
FL2AGM.umap.sym$predicted.id <- as.character(HSC2AGM_CS14.15.query$meta_data$cell_type_pred_knn)
FL2AGM.umap.sym$predicted.score <- as.numeric(HSC2AGM_CS14.15.query$meta_data$cell_type_pred_knn_prob)
FL2AGM.umap.sym <- cbind(FL2AGM.umap.sym, FL_PCW5_HSC.SeuratObj@meta.data)
## Generate plots
cellNum.idx <- table(FL2AGM.umap.sym$predicted.id)
cellNum.idx <- cellNum.idx[c("HSC",
                             "Lympho",
                             "Gr/Prolif",
                             "Mo/Mϕ/Prolif",
                             "Mo/Mϕ")]
### V1 - Plot predicted label
AGM_p3.2_v1 <- ggplot() +
  geom_point(data = bk.umap.AGM, aes(UMAP1, UMAP2, color = Label), 
             inherit.aes = F, alpha = 1, size = 1) +
  geom_point(data = FL2AGM.umap.sym, aes(UMAP1, UMAP2, 
                                         fill = predicted.id), 
             inherit.aes = F, size = 1.2, stroke = .2, shape = 21, alpha = 1) +
  scale_color_manual(values = c("Background" = "#EEEEEE")) +
  scale_fill_manual(values = COL8[names(cellNum.idx)],
                    breaks = names(cellNum.idx),
                    name = "Predicted Labels",
                    labels = paste0(names(cellNum.idx), " (", cellNum.idx, ")")) +
  labs(title = "Project Predicted Labels onto Calvanese et al. data\n(Highlight FL PCW5 HSCs)") +
  guides(color = "none",
         fill = guide_legend(override.aes = list(size = 5))) +
  annotate("segment",
           x = arrow.l$x, xend = arrow.l$x + c(arrow.l$x_len, 0),
           y = arrow.l$y, yend = arrow.l$y + c(0, arrow.l$y_len),
           arrow = arrow(type = "open", length = unit(.1, 'inches')),
           linewidth = .5) +
  annotate("text",
           x = arrow.l$x + .9, y = arrow.l$y + .4, label = "UMAP", 
           size = 3.5,
           color = "black",
           fontface = "bold") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", hjust = .5, colour = "black", size = 13),
        legend.title = element_text(face = "bold", hjust = 0.5, size = 11))
### Save plot
Cairo::CairoPDF("AGM_p3.2_v1.pdf", width = 7, height = 6, family = "Arial")
AGM_p3.2_v1
dev.off()
### V2 - Plot HSC density
nbins = 50
CT.idx = "HSC"
### density version
AGM_p3.2_v2.1 <- ggplot() +
  geom_point(data = bk.umap.AGM, aes(UMAP1, UMAP2, color = Label),
             inherit.aes = F, alpha = 1, size = 1) +
  scale_color_manual(values = c("Background" = "#EEEEEE")) +
  geom_density_2d_filled(data = FL2AGM.umap.sym,
                         aes(x = UMAP1, y = UMAP2, alpha = ..level..),
                         inherit.aes = F, bins = nbins) +
  scale_fill_manual(values = colorRampPalette(c("#F1E7E9", COL8[CT.idx]))(nbins)) +
  scale_alpha_manual(values = c(0, rep(1, nbins-1))) +
  labs(title = "Projection onto Calvanese et al. data\n(Density of FL PCW5 HSCs)") +
  annotate("segment",
           x = arrow.l$x, xend = arrow.l$x + c(arrow.l$x_len, 0),
           y = arrow.l$y, yend = arrow.l$y + c(0, arrow.l$y_len),
           arrow = arrow(type = "open", length = unit(.1, 'inches')),
           linewidth = .5) +
  annotate("text",
           x = arrow.l$x + .9, y = arrow.l$y + .4, label = "UMAP", 
           size = 3.5,
           color = "black",
           fontface = "bold") +
  theme_void() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = .5, colour = "black", size = 13))
### Save plot
Cairo::CairoPDF("AGM_p3.2_v2.1.pdf", width = 6, height = 6, family = "Arial")
AGM_p3.2_v2.1
dev.off()
### contour version
AGM_p3.2_v2.2 <- ggplot() +
  geom_point(data = bk.umap.AGM, aes(UMAP1, UMAP2, color = Label),
             inherit.aes = F, alpha = 1, size = 1) +
  scale_color_manual(values = c("Background" = "#EEEEEE")) +
  stat_density_2d(data = FL2AGM.umap.sym,
                  aes(x = UMAP1, y = UMAP2, fill = ..level.., alpha = ifelse(..level..==0,0,1)),
                  geom = "polygon", contour = TRUE, bins = nbins+2, inherit.aes = F,
                  color = "black", linewidth = 0.07) +
  scale_fill_gradient(low = "#F1E7E9", high = COL8[CT.idx]) +
  labs(title = "Projection Density of FL PCW5 HSCs") +
  annotate("segment",
           x = arrow.l$x, xend = arrow.l$x + c(arrow.l$x_len, 0),
           y = arrow.l$y, yend = arrow.l$y + c(0, arrow.l$y_len),
           arrow = arrow(type = "open", length = unit(.1, 'inches')),
           linewidth = .5) +
  annotate("text",
           x = arrow.l$x + .9, y = arrow.l$y + .4, label = "UMAP", 
           size = 3.5,
           color = "black",
           fontface = "bold") +
  theme_void() +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = .5, colour = "black", size = 13))
### Save plot
Cairo::CairoPDF("AGM_p3.2_v2.2.pdf", width = 6, height = 6, family = "Arial")
AGM_p3.2_v2.2
dev.off()

################
### Bar plot ###
################

# Bar plot of donor v.s cell type
## Preapare plot data
bar_AGM_cross_meta.df <- data.frame(Donor = FL_PCW5_HSC.SeuratObj$donorID, 
                                    Predicted_Labels = as.character(HSC2AGM_CS14.15.query$meta_data$cell_type_pred_knn))
bar_AGM_cross_meta.df %<>%
  group_by(Donor, Predicted_Labels) %>%
  tally() %>%
  mutate(Proportion = n / sum(n))
PL.v <- c("HSC",
          "Lympho",
          "Gr/Prolif",
          "Mo/Mϕ/Prolif",
          "Mo/Mϕ")
bar_AGM_cross_meta.df$Predicted_Labels <- factor(bar_AGM_cross_meta.df$Predicted_Labels, levels = rev(PL.v))
## Generate plot
### Proportion
AGM_p3.3 <- ggplot(bar_AGM_cross_meta.df, aes(x = Donor, y = Proportion, fill = Predicted_Labels)) +
  geom_bar(stat = "identity", position = "stack", width = .5) +
  scale_fill_manual(values = COL8[PL.v],
                    breaks = PL.v,
                    name = "Predicted Labels") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Donor", y = "Proportion", title = "Proportion of Predicted Labels by Donor") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = .5, colour = "black", size = 13),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text = element_text(face = "bold", hjust = .5, colour = "black", size = 9), 
        axis.title = element_text(face = "bold", hjust = .5, colour = "black", size = 12),
        legend.title = element_text(face = "bold", hjust = 0, size = 11))
### Absolute
AGM_p3.3.abs <- ggplot(bar_AGM_cross_meta.df, aes(x = Donor, y = n, fill = Predicted_Labels)) +
  geom_bar(stat = "identity", position = "stack", width = .5) +
  scale_fill_manual(values = COL8[PL.v],
                    breaks = PL.v,
                    name = "Predicted Labels") +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Donor", y = "Cell number", title = "Number of Predicted Labels by Donor") +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = .5, colour = "black", size = 13),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        axis.text = element_text(face = "bold", colour = "black", size = 9), 
        axis.title = element_text(face = "bold", colour = "black", size = 12),
        legend.title = element_text(face = "bold", hjust = 0, size = 11))
AGM_p3.3.all <- cowplot::plot_grid(AGM_p3.3, AGM_p3.3.abs, ncol = 1, labels = NULL)
## Save plot
Cairo::CairoPDF("AGM_p3.3.pdf", width = 9, height = 11, family = "Arial")
AGM_p3.3.all
dev.off()

# sessionInfo
sessionInfo()