#### Set Up ####
library(tidyverse)
library(patchwork)
library(ggplot2)
library(ggplotify)
library(ComplexHeatmap)
library(circlize)
library(openxlsx)
set.seed(777)

workdir <- '/work/Jupyterlab/Project/CellCycle/02.Dynamic_CC/06.GAM_CycSpline/'
datadir <- paste0(workdir,'data/')
plotdir <- paste0(workdir,'plots/')

#### Load Data ####
marker <- readRDS('/work/Jupyterlab/Project/CellCycle/02.Dynamic_CC/utils/CC_marker.rds')
res <- readRDS(paste0(datadir,'10.GEX_fit_CycSpline.rds'))
gene_fit <- res$gene_fit
yhats <- res$yhats

idx <- which(gene_fit$p_adjust < 0.01 & gene_fit$avg_log2fc > quantile(gene_fit$avg_log2fc,0.9))
yhatScaled <- t(scale(t(yhats[idx,])))

gene_fit_sig <- gene_fit[idx,] |> rownames_to_column('Gene')
write.xlsx(gene_fit_sig,paste0(datadir,'11.GEX_gene_fit_sig.xlsx'))

#### K-means clustering ####
clusters <- kmeans(yhatScaled, centers = 5)$cluster
clusterLabels <- clusters
cUniq <- unique(clusterLabels)
cUniq <- sort(cUniq)
p <- list()
nPointsClus <- 40
for (xx in cUniq) {
  clus <- as.character(xx)
  cId <- which(clusterLabels == xx)
  p[[clus]] <- ggplot(data = data.frame(x = 1:nPointsClus,
                                        y = rep(range(yhatScaled[cId, ]),
                                                nPointsClus / 2)),
                      aes(x = x, y = y)) +
    geom_point(alpha = 0) +
    labs(title = paste0("Cluster ", xx),  x = "Pseudotime", y = "Normalized expression") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(breaks = c(0,10,20,30,40), 
                       labels = c('0', expression(pi/2), expression(pi), expression(3*pi/2), expression(2*pi)))
  for (ii in 1:length(cId)) {
    geneId <- rownames(yhatScaled)[cId[ii]]
    p[[clus]] <- p[[clus]] +
      geom_line(data = data.frame(x = rep(1:nPointsClus, 2),
                                  y = yhatScaled[geneId, ],
                                  lineage = rep(0, each = nPointsClus)),
                aes(col = as.character(lineage), group = lineage), lwd = 0.1)
  }
  p[[clus]] <- p[[clus]] + guides(color = FALSE) +
    scale_color_manual(values = c("orange"),
                       breaks = c("0"))  
}

pdf(paste0(plotdir,"11.GEX_PatternCluster_kmeans.pdf"), width=15, height=20)
page1 <- wrap_plots(p[1:5], ncol = 2)
print(page1)#;print(page2);print(page3)
dev.off()

#### Heatmap ####
clusterbyphase <- list()
clusterbyphase[['g1']] <- c(4) 
clusterbyphase[['g1s']] <- c(2) 
clusterbyphase[['s']] <- c(3) 
clusterbyphase[['g2']] <- c(5)
clusterbyphase[['m']] <- c(1)

gex <- lapply(clusterbyphase, function(phase){
  cId <- which(clusterLabels %in% phase)
  gexname <- rownames(yhatScaled[cId,])
  gexname
})
gexs <- unlist(gex)

heatmat = yhatScaled[gexs,]
cc_point <- seq(0,2*pi,length.out = 40)
cc_point <- round(as.numeric(cc_point)/pi,digits = 2)
pi_labels <- ifelse(
  cc_point == 0, "0",
  ifelse(
    cc_point == 1, "pi",
    paste0(cc_point, "*pi")
  )
)

expr_labels <- parse(text = pi_labels)
colnames(heatmat) <- pi_labels

col_fun <- colorRamp2(c(-1.5, 0, 1.5), c("#3a86ff", "white", "#d62828"))
p <- Heatmap(heatmat,
             name = "Z-score",
             col = col_fun,
             heatmap_width  = unit(12,'cm'),
             cluster_rows = FALSE,
             cluster_columns = FALSE,
             show_row_names = FALSE,
             show_column_names = TRUE,
             column_names_gp = grid::gpar(fontsize = 5),
             use_raster = TRUE,
             column_labels = expr_labels,
             heatmap_legend_param = gpar(fontsize = 5, fontface = "bold"),
             raster_quality = 5
)


pt <- p + rowAnnotation(link = anno_mark(at = match(unlist(marker), rownames(heatmat)),
                                         labels = unlist(marker), labels_gp = gpar(fontsize = 5)),width = unit(1, "cm"))

cairo_pdf(file = paste0(plotdir, "11.GEX_CC_heatmap_kmeans.pdf"),
          width = 8, height = 10,family = 'Arial')
draw(pt,
     heatmap_legend_side      = "left",
     annotation_legend_side   = "left",
     padding                  = unit(c(5,20,5,5), "mm"),  # top, right, bottom, left
)
dev.off()

#### order heatmap by gene peak phase ####
order <- gene_fit[idx,] |> arrange(peak_phase)
heatmat = yhatScaled[rownames(order),]

col_fun <- colorRamp2(c(-1.5, 0, 1.5), c("#3a86ff", "white", "#d62828"))
p <- Heatmap(heatmat,
             name = "Z-score",
             col = col_fun,
             heatmap_width  = unit(12,'cm'),
             cluster_rows = FALSE,
             cluster_columns = FALSE,
             show_row_names = FALSE,
             show_column_names = TRUE,
             column_names_gp = grid::gpar(fontsize = 5),
             use_raster = TRUE,
             column_labels = expr_labels,
             raster_quality = 5
)
pt <- p + rowAnnotation(link = anno_mark(at = match(unlist(marker), rownames(heatmat)),
                                         labels = unlist(marker), labels_gp = gpar(fontsize = 5)),width = unit(1, "cm"))

cairo_pdf(file = paste0(plotdir, "11.GEX_CC_heatmap_orderedbypeakphase.pdf"),
          width = 8, height = 10,family = 'Arial')
draw(pt,
     heatmap_legend_side      = "left",
     annotation_legend_side   = "left",
     padding                  = unit(c(5,20,5,5), "mm"),  # top, right, bottom, left
)
dev.off()

# save
saveRDS(gex,paste0(datadir,'11.GEX_clusterbyphase_kmeans.rds'))

#### Add text color ####
gex <- readRDS(paste0(datadir,'11.GEX_clusterbyphase_kmeans.rds'))
gexs <- unlist(gex)
order <- gene_fit[idx,] |> arrange(peak_phase)
heatmat = yhatScaled[rownames(order),]
cc_point <- seq(0,2*pi,length.out = 40)
cc_point <- round(as.numeric(cc_point)/pi,digits = 2)
pi_labels <- ifelse(
  cc_point == 0, "0",
  ifelse(
    cc_point == 1, "pi",
    paste0(cc_point, "*pi")
  )
)

expr_labels <- parse(text = pi_labels)
colnames(heatmat) <- pi_labels

# row annotation
marker$g1<-NULL
gene_df <- data.frame(
  gene = unlist(marker)
) |>
  dplyr::mutate(
    gene_group = case_when(
      gene %in% marker$s   ~ "S",
      gene %in% marker$g2m ~ "G2M",
      TRUE ~ NA_character_
    )
  )
gene_df$gene_group <- factor(gene_df$gene_group, levels = c('S','G2M'))
gene_group_levels <- levels(gene_df$gene_group)
group_colors <- setNames(c('#ff7f0e','#2ca02c'), gene_group_levels)
gene_df$gene_color <- group_colors[as.character(gene_df$gene_group)]  
  
#range(heatmat)
#quantile(heatmat, c(0.01, 0.05, 0.1, 0.5, 0.8, 0.9, 0.95, 0.99), na.rm=TRUE)
col_fun <- colorRamp2(c(-1.5, 0, 1.5), c("#3a86ff", "white", "#d62828"))

ha <- rowAnnotation(foo = anno_mark(at = match(unlist(marker),rownames(heatmat)), 
                                    labels = unlist(marker),
                                    labels_gp = gpar(fontsize = 5,
                                                     col = gene_df$gene_color,
                                                     fontfamily = 'Arial',
                                                     fontface = 'plain'),
                                    link_gp = gpar(col = gene_df$gene_color,lwd=0.5),
                                    link_width = unit(5, "mm"))
)

gene_group_legend <- Legend(
  title = "Phase",
  labels = names(group_colors),
  legend_gp = gpar(fill = group_colors),
  labels_gp = gpar(fontsize = 10),
  title_gp = gpar(fontsize = 10,fontface = "bold"),
  #column_gap = unit(2, "mm"),
  row_gap = unit(2, "mm")
)
pd = packLegend(gene_group_legend, max_width = unit(5.5, "cm"))

p <- Heatmap(heatmat,
             name = "Z-score",
             col = col_fun,
             right_annotation = ha,
             heatmap_width  = unit(13,'cm'),
             cluster_rows = FALSE,
             cluster_columns = FALSE,
             show_row_names = FALSE,
             show_column_names = TRUE,
             column_names_gp = grid::gpar(fontsize = 5),
             use_raster = TRUE,
             column_labels = expr_labels,
             heatmap_legend_param = gpar(fontsize = 5, fontface = "bold"),
             raster_quality = 5
)

cairo_pdf(file = paste0(plotdir, "11.GEX_CC_heatmap_orderedbypeakphase_textcolored.pdf"),
          width = 8, height = 10,family = 'Arial')
draw(p,
     heatmap_legend_side      = "left",
     annotation_legend_side   = "left",
     #padding                  = unit(c(5,10,5,5), "mm"),  # top, right, bottom, left
     annotation_legend_list = list(pd)
)
dev.off()
