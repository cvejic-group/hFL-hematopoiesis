####################################
### Metaplot eRegulon Validation ###
####################################

###-------------------Note-------------------###
###-------      Run with R-4.4.2      -------###
###------------------------------------------###

##################
# Initialization #
##################

setwd("~/")
dir.create(paste0(FIG_DIR, "01.metaplot_Validation_eRegulon"))
dir.create(paste0(RES_DIR, "01.metaplot_Validation_eRegulon"))

# Load addtional packages
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

##################
# Set parameters #
##################

xlsx_path <- "./data/CT_eGRN_DAR_Prune_Metadata.xlsx"
raw_xlsx_path <- "./data/eRegulon_Raw_Metadata.xlsx"
out_dir   <- "./plots/01.metaplot_Validation_eRegulon/"

# Load raw eRegulons
message("Reading raw eRegulon (+/+ only) ...")
raw_df <- read_xlsx(raw_xlsx_path, sheet = "+/+")

# Path to bw files
bw_dir <- "~/RefData/ChIPseq/GSE231422/"

# Metaplots setup
CT <- "HSC"
TF.v <- c("ERG", "FLI1", "GATA2", "LMO2", "LYL1", "RUNX1", "TAL1")

#####################
# Generate metaplot #
#####################

# BigWig naming pattern:
#   ChIP: GSE231422_Coverage_<TF>_HSC.bw
#   IgG:  GSE231422_IgG_all4subfractions.bw
bw_igg    <- file.path(bw_dir, "GSE231422_IgG_all4subfractions.bw")

# Metaplot parameters
extend_bp <- 3000 
bin_bp    <- 50

# Heatmap colour scale
col_fun <- colorRamp2(c(0, 1, 2, 4),
                      c("#3a86ff", "#48cae4", "white", "#d62828"))

# ── Helper: parse "chrN:start-end" Region column → centred GRanges ───────────
parse_regions <- function(TF.df) {
  TF_eRegulon.gr <- TF.df %>%
    mutate(
      chr = str_extract(Region, "^chr[^:]+"),
      start = as.integer(str_extract(Region, "(?<=:)[0-9]+")),
      end = as.integer(str_extract(Region, "(?<=-)[0-9]+"))
    ) %>%
    select(chr, start, end) %>%
    as(., "GRanges") %>%
    sort()
  # Return centre-point GRanges (normalizeToMatrix anchors on this)
  resize(TF_eRegulon.gr, width = 1, fix = "center")
}

# ── Helper: read a BigWig ──
bw_cache <- list()
load_bw <- function(path) {
  if (is.null(bw_cache[[path]])) {
    message("  Loading bigwig: ", basename(path))
    bw_cache[[path]] <<- import(path, format = "BigWig")
  }
  bw_cache[[path]]
}

# ── Pre-load IgG once ─────────────────────────────────────────────────────────
stopifnot(file.exists(bw_igg))
bw_igg_gr <- load_bw(bw_igg)

# ── Main loop: one PDF per TF ─────────────────────────────────────────────────
for (tf in TF.v) {
  
  message("\n── Processing TF: ", tf, " ──────────────────────────────────────")
  
  # 1. BigWig path for this TF
  bw_chip_path <- file.path(bw_dir,
                            paste0("GSE231422_Coverage_", tf, "_HSC.bw"))
  
  if (!file.exists(bw_chip_path)) {
    warning("BigWig not found, skipping: ", bw_chip_path)
    next
  }
  
  # 2. Subset eRegulon regions for this TF
  tf_rows <- raw_df %>%
    filter(TF == tf)
  if (nrow(tf_rows) == 0) {
    warning("No eRegulon regions found for TF: ", tf, " in sheet: ", CT)
    next
  }
  message("  eRegulon regions: ", nrow(tf_rows))
  
  # Remove duplicate regions (same region can link to multiple genes)
  unique_regions <- unique(tf_rows$Region)
  message("  Unique regions: ", length(unique_regions))
  
  # 3. Parse regions to centred GRanges
  target_gr <- parse_regions(tf_rows)
  
  # 4. Load ChIP bigwig and compute signal matrices
  bw_chip_gr <- load_bw(bw_chip_path)
  
  message("  Computing ChIP signal matrix ...")
  mat_chip <- normalizeToMatrix(
    signal       = bw_chip_gr,
    target       = target_gr,
    extend       = extend_bp,
    w            = bin_bp,
    value_column = "score",
    background   = 0,
    smooth       = TRUE,
    verbose      = FALSE
  )
  
  message("  Computing IgG signal matrix ...")
  mat_igg <- normalizeToMatrix(
    signal       = bw_igg_gr,
    target       = target_gr,
    extend       = extend_bp,
    w            = bin_bp,
    value_column = "score",
    background   = 0,
    smooth       = TRUE,
    verbose      = FALSE
  )
  
  ## Compute shared y-axis limits from ChIP signal (so IgG looks flat by comparison)
  chip_mean <- colMeans(mat_chip, na.rm = TRUE)
  y_max <- max(chip_mean) * 1.1   # 5% headroom
  # y_min <- min(chip_mean) * 0.95
  shared_ylim <- c(1, y_max)
  
  # 5. Build EnrichedHeatmap objects
  #    Row order: sort by ChIP signal in central ±500 bp
  box_h <- 4
  heat_h <- 9
  ht_chip <- EnrichedHeatmap(
    mat_chip,
    name            = "Signal",
    col             = col_fun,
    column_title    = tf,
    column_title_gp = gpar(fontsize = 13, fontface = "bold"),
    row_title       = paste0(tf, " eRegulon regions\n(n = ", length(unique_regions), ")"),
    row_title_gp    = gpar(fontsize = 10),
    row_order       = order(rowMeans(mat_chip[, (ncol(mat_chip)/2 - 10):(ncol(mat_chip)/2 + 10)]),
                            decreasing = TRUE),
    top_annotation  = HeatmapAnnotation(
      enriched = anno_enriched(
        gp             = gpar(col = "#1F3A6E", lwd = 2),
        axis_param     = list(side = "left", facing = "outside",
                              gp = gpar(fontsize = 8)),
        ylim       = shared_ylim, 
      ),
      height = unit(box_h, "cm")
    ),
    axis_name       = c(paste0("-", extend_bp/1000, " kb"), "0",
                        paste0("+", extend_bp/1000, " kb")),
    axis_name_gp    = gpar(fontsize = 9),
    heatmap_legend_param = list(
      title          = "Intensity",
      title_gp       = gpar(fontsize = 9),
      labels_gp      = gpar(fontsize = 8)
    ),
    height = unit(heat_h, "cm"),
    use_raster      = TRUE,
    raster_quality  = 3
  )
  
  ht_igg <- EnrichedHeatmap(
    mat_igg,
    name                 = "IgG",
    col                  = col_fun,
    column_title         = "IgG",
    column_title_gp      = gpar(fontsize = 13, fontface = "bold"),
    top_annotation       = HeatmapAnnotation(
      enriched = anno_enriched(
        gp = gpar(col = "#1F3A6E", lwd = 2),
        axis = FALSE,
        ylim = shared_ylim, 
      ),
      height = unit(box_h, "cm")
    ),
    axis_name     = c(paste0("-", extend_bp/1000, " kb"), "0",
                      paste0("+", extend_bp/1000, " kb")),
    height = unit(heat_h, "cm"),
    axis_name_gp  = gpar(fontsize = 9),
    show_heatmap_legend = FALSE,
    use_raster    = TRUE,
    raster_quality = 3
  )
  
  # 6. Draw and save
  out_pdf <- paste0(out_dir, paste0(tf, "_HSC_chipseq_metaplot.pdf"))
  pdf(out_pdf, width = 5.5, height = 6.3)
  
  draw(
    ht_chip + ht_igg,
    ht_gap = unit(5, "mm"),
    padding = unit(c(7, 5, 10, 5), "mm")
  )
  grid.text(
    "ChIPseq validation",
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

message("\nDone. All metaplots written to: ", out_dir)
