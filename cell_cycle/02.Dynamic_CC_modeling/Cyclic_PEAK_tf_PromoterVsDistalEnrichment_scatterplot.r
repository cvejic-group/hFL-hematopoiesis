# ============================================================
# Title: Cyclic_PEAK_tf_PromoterVsDistalEnrichment_scatterplot.r
# Purpose:
#   Categorize TF enriched into promoter- and distal-specific
#   Plot scatter plot colored by categories 

# ============================================================
# Load packages
# ============================================================
library(dplyr)
library(stringr)
library(ggplot2)
library(GenomicRanges)
library(SummarizedExperiment)
library(tidyverse)

workdir <- '/work/Jupyterlab/Project/CellCycle/02.Dynamic_CC/06.GAM_CycSpline/'
datadir <- paste0(workdir,'data/')
plotdir <- paste0(workdir,'plots/')

# ============================================================
# Load monalisa sce
# ============================================================
se_promoter <- readRDS(paste0(datadir,'24.Peak_promoterwithin1kb_monalisa.rds'))
se_distal <- readRDS(paste0(datadir,'24.Peak_distalintergenic_monalisa.rds'))
se_all <- readRDS(paste0(datadir,'24.Peak_clusteredbyphase_monalisa.rds'))

# ============================================================
# Categorize and plot
# ============================================================
se_all1 <- apply(assay(se_all, "negLog10Padj"), 1, 
                 function(x) max(abs(x), 0, na.rm = TRUE)) > 3
motif2keep <- names(se_all1[se_all1==TRUE])


sel_promoter <- apply(assay(se_promoter, "negLog10Padj"), 1, 
                      function(x) max(abs(x), 0, na.rm = TRUE)) > 3
sel_distal <- apply(assay(se_distal, "negLog10Padj"), 1, 
                    function(x) max(abs(x), 0, na.rm = TRUE)) > 3
sum(sel_promoter);sum(sel_distal)
motifid2p <- union(names(sel_promoter[sel_promoter==TRUE]),
                   names(sel_distal[sel_distal==TRUE]))
motifname <- as.data.frame(rowData(se_promoter)[motifid2p,c('motif.id','motif.name')])

pro_enr <- se_promoter@assays@data$log2enr[motifid2p,]
pro_pva <- se_promoter@assays@data$negLog10Padj[motifid2p,]
dis_enr <- se_distal@assays@data$log2enr[motifid2p,]
dis_pva <- se_distal@assays@data$negLog10Padj[motifid2p,]

enr_cutoff <- 0
sig_cutoff <- 3
d4p <- list()
d4p_ <- list()
#filter_mode    <- "any" 
for(phase in colnames(pro_enr)){
  phase_pro_enr <- pro_enr |> as.data.frame() |> rownames_to_column('motif.id') |> 
    dplyr::mutate(enr_pro = !!sym(phase)) |>
    dplyr::select(motif.id, enr_pro)
  phase_dis_enr <- dis_enr |> as.data.frame() |> rownames_to_column('motif.id') |> 
    dplyr::mutate(enr_dis = !!sym(phase)) |>
    dplyr::select(motif.id, enr_dis)
  phase_enr <- phase_pro_enr |> left_join(phase_dis_enr, by = 'motif.id')
  
  phase_pro_pva <- pro_pva |> as.data.frame() |> rownames_to_column('motif.id') |> 
    dplyr::mutate(negLog10P_pro = !!sym(phase)) |>
    dplyr::select(motif.id, negLog10P_pro)
  phase_dis_pva <- dis_pva |> as.data.frame() |> rownames_to_column('motif.id') |> 
    dplyr::mutate(negLog10P_dis = !!sym(phase)) |>
    dplyr::select(motif.id, negLog10P_dis)
  phase_pva <- phase_pro_pva |> left_join(phase_dis_pva, by = 'motif.id')
  
  phase_d4p <- phase_enr |> left_join(phase_pva, by = 'motif.id') |>
    left_join(motifname,by='motif.id')
  phase_d4p$nlp_mean <- (phase_d4p$negLog10P_pro + phase_d4p$negLog10P_dis) / 2
  
  # significance per region; direction from sign(log2enr)
  sig_prom <- phase_d4p$negLog10P_pro >= sig_cutoff
  sig_dist <- phase_d4p$negLog10P_dis >= sig_cutoff
  up_prom  <- phase_d4p$enr_pro > enr_cutoff
  up_dist  <- phase_d4p$enr_dis > enr_cutoff
  
  phase_d4p$class <- dplyr::case_when(
    sig_prom &  sig_dist &  up_prom &  up_dist        ~ "Enriched in both for phase",
    sig_prom &  sig_dist & !up_prom & !up_dist        ~ "Enriched in both for other phases",
    sig_prom &  sig_dist & (up_prom != up_dist)       ~ "Opposite",
    sig_prom & !sig_dist                              ~ "Promoter-specific",
    !sig_prom &  sig_dist                             ~ "Distal-specific",
    TRUE                                              ~ "Other"
  )
  
  phase_d4p$in_all <- phase_d4p$motif.id %in% motif2keep
  #phase_d4p$kept <- sig_prom | sig_dist
  class_levels <- c("Enriched in both for phase",
                    "Enriched in both for other phases",
                    "Promoter-specific",
                    "Distal-specific",
                    "Opposite",
                    "Other")
  
  phase_d4p$phase <- phase
  phase_d4p_plot <- phase_d4p
  phase_d4p_plot$class <- factor(phase_d4p_plot$class, levels = class_levels)
  phase_d4p_plot <- phase_d4p_plot[order(phase_d4p_plot$in_all),]
  print(dim(phase_d4p_plot))
  d4p[[phase]] <- phase_d4p
  d4p_[[phase]] <- phase_d4p_plot
}

d4ps <- do.call(rbind, d4p)
d4ps$phase <- factor(d4ps$phase,levels = c('mg1','g1','g1s','s','sg2','g2m'))

lim <- max(abs(c(d4ps$enr_pro, d4ps$enr_dis)), na.rm = TRUE) * 1.1

class_colors <- c("Enriched in both for phase"        = "#E41A1C",
                  "Enriched in both for other phases" = "#377EB8",
                  "Promoter-specific"                 = "#FF7F00",
                  "Distal-specific"                   = "#E41A1C",
                  "Opposite"                          = "#4DAF4A")

p <- ggplot(d4ps[!is.na(d4ps$class), ], aes(enr_pro, enr_dis)) +
  facet_wrap(~phase) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey40") +
  geom_point(aes(size = nlp_mean, fill = class),
             shape = 21, color = "black", stroke = 0.3, alpha = 0.9) +
  ggrepel::geom_text_repel(
    data = subset(d4ps, in_all & class!= 'Other'),
    aes(label = motif.name, color = class),
    size = 2.9, fontface = "bold",
    max.overlaps = Inf, min.segment.length = 0,
    segment.color = "grey70",
    box.padding = 0.5, point.padding = 0.5, force = 0.5,
    show.legend = FALSE
  ) +
  scale_size_continuous(name = "mean -log10 adj.p", range = c(1.5, 6.5)) +
  scale_fill_manual (values = class_colors, name = "Class", drop = TRUE) +
  scale_color_manual(values = class_colors, guide = "none") +
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(x = "log2 enrichment (Up vs Dn) -- Promoter (<=1kb)",
       y = "log2 enrichment (Up vs Dn) -- Distal Intergenic",
       title    = "TF motif enrichment: Promoter vs Distal",
       # subtitle = sprintf("%d motifs (-log10 adj.p > %g at least in one phase)",
       #                    length(unique(d4ps$motif.id)), sig_cutoff)
  ) +
  guides(fill = guide_legend(override.aes = list(size = 5),
                             title.position = "top", ncol = 1),
         size = guide_legend(title.position = "top", ncol = 1)) +
  theme_classic(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold", hjust = 0, size = 13),
        plot.subtitle    = element_text(face = "bold", hjust = 0, size = 9),
        legend.position  = "bottom",
        legend.box       = "horizontal")

cairo_pdf(paste0(plotdir,'28.Peak_clusteredbyphase_monalisa_provsdis_scatter_v3.pdf'),
          width = 20,height = 14,family = 'Arial')
p
dev.off()

# ============================================================
# Write table
# ============================================================
write_motif_enrichment_table <- function(sce){
  motifid2name <- rowData(sce) |> as.data.frame() |>
    select(motif.id, motif.name) 
  
  negLog10Padj <- assay(sce,"negLog10Padj")
  colnames(negLog10Padj) <- paste0(colnames(negLog10Padj),"_negLog10Padj")
  negLog10Padj <- negLog10Padj |> as.data.frame() |> rownames_to_column('motif.id')
  
  log2enr <- assay(sce,"log2enr")
  colnames(log2enr) <- paste0(colnames(log2enr),"_log2enr")
  log2enr <- log2enr |> as.data.frame() |> rownames_to_column('motif.id') 
  
  df2e <- negLog10Padj |>
    left_join(log2enr, by = "motif.id") |>
    left_join(motifid2name, by = "motif.id") |>
    tibble::column_to_rownames("motif.id") |>
    mutate(
      negLog10Padj_sum = rowSums(
        pick(all_of(colnames(negLog10Padj)[-1])),
        na.rm = TRUE
      )
    ) |>
    arrange(desc(negLog10Padj_sum)) |>
    select(-negLog10Padj_sum)
}

promoter_tab <- write_motif_enrichment_table(se_promoter)
distal_tab <- write_motif_enrichment_table(se_distal)

write.csv(promoter_tab, paste0(datadir,'28.Peak_separtedbyPromoter&Distal_scattertable_promoter.csv'))
write.csv(distal_tab, paste0(datadir,'28.Peak_separtedbyPromoter&Distal_scattertable_distal.csv'))
