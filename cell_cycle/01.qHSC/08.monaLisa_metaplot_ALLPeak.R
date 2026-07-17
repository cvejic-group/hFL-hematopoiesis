# ============================================================================
#  ALL Peak metaplot to show accessiablity
# ============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(TFBSTools)
  library(monaLisa)
  library(SummarizedExperiment)
  library(ComplexHeatmap)
  library(BiocParallel)
  library(dplyr)
})

# ---- 0. USER PARAMETERS ----------------------------------------------------
setwd("~/")
res_dir <- "./results/08.monaLisa_metaplot/HSC_all_peaks/"
fig_dir <- "./plots/08.monaLisa_metaplot/HSC_all_peaks/"
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
txdb   <- TxDb.Hsapiens.UCSC.hg38.knownGene
genome <- BSgenome.Hsapiens.UCSC.hg38
species_tax_id <- 9606

# ---- 1. READ HSC PEAKS -----------------------------------------------------
library(data.table)

bed <- fread(
  "~/HSC.bed",
  skip   = 1,       # skip header line
  header = FALSE
)

HSC.gr <- GRanges(
  seqnames = bed$V1,
  ranges   = IRanges(start = bed$V2, end = bed$V3)
)
message("Total HSC peaks: ", length(HSC.gr))

# ---- 2. CHIPSEEKER ANNOTATION ----------------------------------------------
anno_HSC <- annotatePeak(HSC.gr, tssRegion = c(-3000, 3000),
                         TxDb = txdb, level = "transcript",
                         verbose = FALSE)
gr_HSC <- as.GRanges(anno_HSC)

# Subset by annotation
region_mask  <- function(gr, label) startsWith(gr$annotation, label)

HSC_all.gr      <- gr_HSC                                           # all peaks
HSC_promoter.gr <- gr_HSC[region_mask(gr_HSC, "Promoter (<=1kb)")] # promoter
HSC_distal.gr   <- gr_HSC[region_mask(gr_HSC, "Distal Intergenic")]# distal

message("All: ",      length(HSC_all.gr),
        "  Promoter: ", length(HSC_promoter.gr),
        "  Distal: ",   length(HSC_distal.gr))

saveRDS(list(HSC_all.gr = HSC_all.gr,
             HSC_promoter.gr = HSC_promoter.gr,
             HSC_distal.gr = HSC_distal.gr), 
        file = paste0(res_dir, "HSC_peak_annotation.rds"))

# ---- 3. Prepare metaplot v1 ---------------------------------------------
library(caret)
library(pROC)
library(future)
library(furrr)
library(readxl)
library(soGGi)
library(ggrepel)

# If only BAM is available, convert first
qHSC_BAM <- "~/HSC-Q/HSC-Q.dedup.bam"
cHSC_BAM <- "~/HSC-C/HSC-C.dedup.bam"
qHSC_bw <- "~/HSC.Q-TileSize-50-normMethod-ReadsInTSS-ArchR.bw"
cHSC_bw <- "~/HSC.C-TileSize-50-normMethod-ReadsInTSS-ArchR.bw"

## Compute profile for HSC — three region sets
### All HSC peaks
HSC_all.list <- list()
HSC_all.list[[1]] <- regionPlot(qHSC_bw, HSC_all.gr,
                                style = "point", format = "bigwig")
HSC_all.list[[2]] <- regionPlot(cHSC_bw, HSC_all.gr,
                                style = "point", format = "bigwig")
names(HSC_all.list) <- c("qHSC", "non-qHSC")

### Promoter peaks
HSC_promoter.list <- list()
HSC_promoter.list[[1]] <- regionPlot(qHSC_bw, HSC_promoter.gr,
                                     style = "point", format = "bigwig")
HSC_promoter.list[[2]] <- regionPlot(cHSC_bw, HSC_promoter.gr,
                                     style = "point", format = "bigwig")
names(HSC_promoter.list) <- c("qHSC", "non-qHSC")

### Distal peaks
HSC_distal.list <- list()
HSC_distal.list[[1]] <- regionPlot(qHSC_bw, HSC_distal.gr,
                                   style = "point", format = "bigwig")
HSC_distal.list[[2]] <- regionPlot(cHSC_bw, HSC_distal.gr,
                                   style = "point", format = "bigwig")
names(HSC_distal.list) <- c("qHSC", "non-qHSC")

# ---- 4. Plot all four conditions in a loop ---------------------------------

plot_configs <- list(
  list(data    = HSC_all.list,
       title   = "Coverage Plot: All HSC Peaks (qHSC vs. non-qHSC)",
       outfile = "HSC_metaplot_All.pdf"),
  list(data    = HSC_promoter.list,
       title   = "Coverage Plot: HSC Promoter Peaks (qHSC vs. non-qHSC)",
       outfile = "HSC_metaplot_Promoter.pdf"),
  list(data    = HSC_distal.list,
       title   = "Coverage Plot: HSC Distal Peaks (qHSC vs. non-qHSC)",
       outfile = "HSC_metaplot_Distal.pdf")
)

COL2 <- c("#c1121f", "#669bbc")
CT.v <- c("qHSC", "non-qHSC")

for (cfg in plot_configs) {
  
  # Extract signal matrix for each cell type
  signal_list <- lapply(CT.v, function(ct) {
    mat    <- assays(cfg$data[[ct]])[[1]]
    extend <- (ncol(mat) - 1) / 2
    pos    <- seq(-extend, extend, length.out = ncol(mat))
    
    data.frame(
      position = pos,
      coverage = colMeans(mat, na.rm = TRUE),   # BigWig already normalised by ArchR
      CT       = ct
    )
  })
  
  signal_df        <- bind_rows(signal_list)
  signal_df$CT     <- factor(signal_df$CT, levels = c("qHSC", "non-qHSC"))
  
  label_df <- signal_df %>%
    group_by(CT) %>%
    filter(position >= -200 & position <= 200) %>%
    slice_max(coverage, n = 1) %>%
    ungroup()
  
  # Dynamic nudge based on actual data range
  x_range  <- diff(range(signal_df$position, na.rm = TRUE))
  y_range  <- diff(range(signal_df$coverage, na.rm = TRUE))
  nudge_x  <- x_range * 0.15    # 15% of x range
  nudge_y  <- y_range * 0.1     # 10% of y range
  
  p <- ggplot(signal_df, aes(x = position, y = coverage, colour = CT)) +
    geom_line(linewidth = 1.5, na.rm = TRUE) +
    geom_text_repel(
      data           = label_df,
      aes(label      = CT),
      size           = 3,
      fontface       = "bold",
      direction      = "y",
      nudge_x        = nudge_x,
      nudge_y        = nudge_y,
      segment.size   = 0.3,
      segment.colour = "grey50",
      show.legend    = FALSE,
      max.overlaps   = Inf
    ) +
    scale_color_manual(values = COL2) +
    labs(
      x      = "Distance from center (bp)",
      y      = "Mean accessibility (ReadsInTSS normalised)",
      title  = cfg$title,
      colour = NULL
    ) +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", hjust = .5,
                                    colour = "black", size = 13))
  
  Cairo::CairoPDF(paste0(fig_dir, cfg$outfile),
                  width = 8, height = 6, family = "Arial")
  print(p)
  dev.off()
  message("Saved: ", cfg$outfile)
}

# ---- 5. Generate metaplot V2 ---------------------------------------------
library(Rsamtools)
library(GenomicAlignments)
library(rtracklayer)
library(caret)
library(pROC)
library(future)
library(furrr)
library(readxl)
library(soGGi)
library(EnrichedHeatmap)
library(circlize)
library(ComplexHeatmap)
library(rtracklayer)

atac_configs <- list(
  list(gr      = HSC_all.gr,
       label   = "HSC All Peaks",
       outfile = "HSC_All_metaplot_v2.pdf"),
  list(gr      = HSC_promoter.gr,
       label   = "HSC Promoter Peaks",
       outfile = "HSC_Promoter_metaplot_v2.pdf"),
  list(gr      = HSC_distal.gr,
       label   = "HSC Distal Peaks",
       outfile = "HSC_Distal_metaplot_v2.pdf")
)

# Load BigWigs once
bw_qHSC_gr <- import(qHSC_bw, format = "BigWig")
bw_cHSC_gr <- import(cHSC_bw, format = "BigWig")

# Shared parameters
extend_bp <- 3000
bin_bp    <- 50
box_h     <- 4
heat_h    <- 9
col_fun_atac <- colorRamp2(c(0, 0.5, 1, 1.5, 2),
                           # c("#3a86ff", "#48cae4", "white", "#d62828")
                           c("#F9F6F4", "#F1DEDA", "#C684AA", "#6B4E80", "#4D4665")
)

# ── Main loop ─────────────────────────────────────────────────────────────────
for (cfg in atac_configs) {
  
  message("\n── Processing: ", cfg$label, " ──")
  
  # Centre the target regions
  target_gr <- resize(cfg$gr, width = 1, fix = "center")
  n_regions <- length(target_gr)
  
  # Compute signal matrices
  message("  normalizeToMatrix: qHSC ...")
  mat_qHSC <- normalizeToMatrix(
    signal       = bw_qHSC_gr,
    target       = target_gr,
    extend       = extend_bp,
    w            = bin_bp,
    value_column = "score",
    background   = 0,
    smooth       = TRUE,
    verbose      = FALSE
  )
  
  message("  normalizeToMatrix: non-qHSC ...")
  mat_cHSC <- normalizeToMatrix(
    signal       = bw_cHSC_gr,
    target       = target_gr,
    extend       = extend_bp,
    w            = bin_bp,
    value_column = "score",
    background   = 0,
    smooth       = TRUE,
    verbose      = FALSE
  )
  
  # Shared y-axis from qHSC signal
  qHSC_mean  <- colMeans(mat_qHSC, na.rm = TRUE)
  cHSC_mean  <- colMeans(mat_cHSC, na.rm = TRUE)
  shared_ylim <- c(0, max(qHSC_mean, cHSC_mean) * 1.1)
  
  # Dynamic colour scale based on both qHSC and cHSC signal range
  signal_max  <- max(
    quantile(mat_qHSC, 0.99, na.rm = TRUE),
    quantile(mat_cHSC, 0.99, na.rm = TRUE)
  )
  signal_mid  <- signal_max * 0.5
  
  col_fun_atac <- colorRamp2(
    c(0, signal_mid, signal_max),
    c("#3a86ff", "white", "#d62828")
  )
  
  # Row order by central signal
  central_cols <- (ncol(mat_qHSC)/2 - 10):(ncol(mat_qHSC)/2 + 10)
  row_ord      <- order(rowMeans(mat_qHSC[, central_cols]), decreasing = TRUE)
  
  # ── qHSC heatmap ──────────────────────────────────────────────────────────
  ht_qHSC <- EnrichedHeatmap(
    mat_qHSC,
    name            = "Accessibility",
    col             = col_fun_atac,
    column_title    = "qHSC",
    column_title_gp = gpar(fontsize = 13, fontface = "bold"),
    row_title       = paste0(cfg$label, "\n(n = ", n_regions, ")"),
    row_title_gp    = gpar(fontsize = 10),
    row_order       = row_ord,
    top_annotation  = HeatmapAnnotation(
      enriched = anno_enriched(
        gp         = gpar(col = "#1F3A6E", lwd = 2),
        ylim       = shared_ylim,
        axis_param = list(side = "left", facing = "outside",
                          gp = gpar(fontsize = 8))
      ),
      height = unit(box_h, "cm")
    ),
    axis_name       = c(paste0("-", extend_bp/1000, " kb"), "0",
                        paste0("+", extend_bp/1000, " kb")),
    axis_name_gp    = gpar(fontsize = 9),
    heatmap_legend_param = list(
      title     = "Accessibility",
      title_gp  = gpar(fontsize = 9),
      labels_gp = gpar(fontsize = 8)
    ),
    height          = unit(heat_h, "cm"),
    use_raster      = TRUE,
    raster_quality  = 3
  )
  
  # ── non-qHSC heatmap ──────────────────────────────────────────────────────
  ht_cHSC <- EnrichedHeatmap(
    mat_cHSC,
    name                 = "non-qHSC",
    col                  = col_fun_atac,
    column_title         = "non-qHSC",
    column_title_gp      = gpar(fontsize = 13, fontface = "bold"),
    row_order            = row_ord,
    top_annotation       = HeatmapAnnotation(
      enriched = anno_enriched(
        gp         = gpar(col = "#1F3A6E", lwd = 2),
        ylim       = shared_ylim,
        axis       = FALSE
      ),
      height = unit(box_h, "cm")
    ),
    axis_name     = c(paste0("-", extend_bp/1000, " kb"), "0",
                      paste0("+", extend_bp/1000, " kb")),
    axis_name_gp  = gpar(fontsize = 9),
    height        = unit(heat_h, "cm"),
    show_heatmap_legend = FALSE,
    use_raster    = TRUE,
    raster_quality = 3
  )
  
  # ── Draw and save ──────────────────────────────────────────────────────────
  out_pdf <- paste0(fig_dir, cfg$outfile)
  pdf(out_pdf, width = 5.5, height = 6.3)
  
  draw(
    ht_qHSC + ht_cHSC,
    ht_gap  = unit(5, "mm"),
    padding = unit(c(7, 5, 10, 5), "mm")
  )
  
  grid.text(
    "scATAC-seq Accessibility",
    x    = 0.5,
    y    = unit(1, "npc") - unit(3, "mm"),
    just = c("centre", "top"),
    gp   = gpar(fontsize = 15, fontface = "bold")
  )
  
  grid.text(
    "Distance to Region Center",
    x    = 0.5,
    y    = unit(5, "mm"),
    just = c("centre", "bottom"),
    gp   = gpar(fontsize = 11)
  )
  
  dev.off()
  message("  Saved: ", out_pdf)
}

