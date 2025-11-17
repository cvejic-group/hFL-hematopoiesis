#### Set Up ####
library(tidyverse)
library(patchwork)
library(ggplot2)
library(ggplotify)
library(ComplexHeatmap)
library(circlize)
set.seed(777)

workdir <- '/work/Jupyterlab/Project/CellCycle/02.Dynamic_CC/06.GAM_CycSpline/'
datadir <- paste0(workdir,'data/')
plotdir <- paste0(workdir,'plots/')

#### Load Data ####
res <- readRDS(paste0(datadir,'20.Peak_fit_CycSpline.rds'))
gene_fit <- res$gene_fit
yhats <- res$yhats

idx <- which(gene_fit$p_adjust < 0.01 & gene_fit$avg_log2fc > quantile(gene_fit$avg_log2fc,0.85))
yhatScaled <- t(scale(t(yhats[idx,])))

peak_fit_sig <- gene_fit[idx,] |> rownames_to_column('Peak')
write.xlsx(peak_fit_sig,paste0(datadir,'21.Peak_gene_fit_sig.xlsx'))

#### order heatmap by gene peak phase ####
order <- gene_fit[idx,] |> arrange(peak_phase)
heatmat = yhatScaled[rownames(order),]
timepoint <- seq(0,2*pi,length.out=40)
timepoint <- round(timepoint/pi,digits = 2)
timepoint <- ifelse(
  timepoint == 0, "0",
  ifelse(
    timepoint ==1 ,'pi',
    paste0(timepoint,'*pi')
  )
)
expr_labels <- parse(text = timepoint)
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


cairo_pdf(file = paste0(plotdir, "21.Peak_CC_heatmap_orderedbypeakphase.pdf"),
          width = 8, height = 10,family = 'Arial')
draw(p,
     heatmap_legend_side      = "left",
     annotation_legend_side   = "left",
     padding                  = unit(c(5,20,5,5), "mm"),  # top, right, bottom, left
)
dev.off()

#### K-means clustering ####
clusters <- kmeans(yhatScaled, centers = 12)$cluster
clusterLabels <- clusters
cUniq <- unique(clusterLabels)
cUniq <- sort(cUniq)
p <- list()
nPointsClus <- 40

for (xx in cUniq) {
  clus <- as.character(xx)
  cId  <- which(clusterLabels == xx)
  
  p[[clus]] <- ggplot(
    data = data.frame(
      x = 1:nPointsClus,
      y = rep(range(yhatScaled[cId, ], na.rm = TRUE), nPointsClus / 2)
    ),
    aes(x = x, y = y)
  ) +
    geom_point(alpha = 0) +
    labs(title = paste0("Cluster ", xx),
         x = "Pseudotime", y = "Normalized expression") +
    theme(plot.title = element_text(hjust = 0.5),
          panel.background = element_blank(),
          panel.border = element_rect(fill=NA,color='black'),
          panel.grid = element_line(color='lightgrey',linetype='dotted')) +
    scale_x_continuous(
      breaks = c(0, 10, 20, 30, 40),
      labels = c('0', expression(pi/2), expression(pi), expression(3*pi/2), expression(2*pi))
    )
  
  mat <- yhatScaled[cId, , drop = FALSE]
  df_lines <- data.frame(
    x      = rep(1:nPointsClus, times = length(cId)),
    y      = as.vector(t(mat)),
    gene   = rep(rownames(yhatScaled)[cId], each = nPointsClus),
    lineage = 0
  )
  
  p[[clus]] <- p[[clus]] +
    geom_line(data = df_lines,
              aes(x = x, y = y, group = gene, col = as.character(lineage)),
              linewidth = 0.1) +
    guides(color = FALSE) +
    scale_color_manual(values = c("orange"), breaks = c("0"))
}


pdf(paste0(plotdir,"21.Peak_PatternCluster_kmeans.pdf"), width=20, height=25)
page1 <- wrap_plots(p[1:12], ncol = 3)
print(page1)#;print(page2);print(page3)
dev.off()

#### Heatmap ####
clusterbyphase <- list()
clusterbyphase[['mg1']] <- c(11) 
clusterbyphase[['g1']] <- c(8,9) 
clusterbyphase[['g1s']] <- c(6,12,10) 
clusterbyphase[['s']] <- c(4,7) 
clusterbyphase[['sg2']] <- c(1,2) 
clusterbyphase[['g2m']] <- c(5,3)



peak <- lapply(clusterbyphase, function(phase){
  cId <- which(clusterLabels %in% phase)
  gexname <- rownames(yhatScaled[cId,])
  gexname
})
peaks <- unlist(peak)

heatmat = yhatScaled[peaks,]
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



cairo_pdf(file = paste0(plotdir, "21.Peak_CC_heatmap_kmeans.pdf"),
          width = 8, height = 10,family = 'Arial')
draw(p,
     heatmap_legend_side      = "left",
     annotation_legend_side   = "left",
     padding                  = unit(c(5,20,5,5), "mm"),  # top, right, bottom, left
)
dev.off()


# save
saveRDS(peak, paste0(datadir,'21.Peak_clusterbyphase_kmeans.rds'))

