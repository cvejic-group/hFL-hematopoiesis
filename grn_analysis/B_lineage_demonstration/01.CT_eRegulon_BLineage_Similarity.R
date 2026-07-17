#########################################
### CT-eRegulon B Lineage Similarity  ###
#########################################

###-------------------Note-------------------###
###-------      Run with R-4.4.2      -------###
###------------------------------------------###

##################
# Initialization #
##################

setwd("~/")
source("./00.Initialization.R")
dir.create(paste0(FIG_DIR, "01.CT_eRegulon_BLineage_Similarity"))

# Load addtional packages
library(caret)
library(pROC)
library(future)
library(furrr)
library(readxl)

#############################################
# Compute B Lineage CT-eRegulon similarity #
#############################################

xlsx_path <- "./data/CT_eGRN_DAR_Prune_Metadata.xlsx"   # per-cell-type eRegulon metadata

# TFs to compare
TF.v <- c("EBF1", "PAX5", "TCF3", "FOXO1")

# Cell types to plot
CT.v <- c("HSC", "PreProB", 
          # "ProB-1", 
          "ProB-2",
          # "Large-PreB",
          "Small-PreB", 
          "IM-B")

# Output file
out_pdf <- "./plots/01.CT_eRegulon_BLineage_Similarity/JI_heatmap.pdf"
out_pdf_norm <- "./plots/01.CT_eRegulon_BLineage_Similarity/JI_heatmap_normalised.pdf"

# Figure dimensions (inches): 2 columns (gene | region), rows = number of CTs
fig_width  <- 10
fig_height <- 5 * length(CT.v)
# ══════════════════════════════════════════════════════════════════════════════

# ── Helper: compute pairwise Jaccard Index matrix ────────────────────────────
#
# @param df     data.frame for one cell type
# @param tf_vec character vector of TF names; missing TFs get JI = 0
# @param col    column name for set elements ("Gene" or "Region")
#
# @return data.frame in long format: TF_x, TF_y, JI
compute_JI <- function(df, tf_vec, col = "Gene") {
  
  # Build a named list: TF -> unique set of targets (empty set if TF absent)
  sets <- lapply(tf_vec, function(tf) unique(df[[col]][df$TF == tf]))
  names(sets) <- tf_vec
  
  # Compute pairwise JI; missing TF pairs receive 0
  n <- length(tf_vec)
  ji_mat <- matrix(NA_real_, n, n, dimnames = list(tf_vec, tf_vec))
  
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      A <- sets[[i]]; B <- sets[[j]]
      inter <- length(intersect(A, B))
      union <- length(union(A, B))
      ji_mat[i, j] <- if (union == 0) 0 else inter / union
    }
  }
  
  as.data.frame(as.table(ji_mat), stringsAsFactors = FALSE) |>
    setNames(c("TF_x", "TF_y", "JI"))
}


# ── Helper: build one heatmap ggplot object ───────────────────────────────────
#
# @param ji_long     long-format data.frame with columns TF_x, TF_y, and value_col
# @param tf_vec      ordered TF names (sets axis order)
# @param title       string, panel title
# @param value_col   column to use for fill ("JI" for raw, "JI_norm" for normalised)
# @param legend_name legend title string
# @param show_legend logical
#
# @return ggplot object
make_heatmap <- function(ji_long, tf_vec, title,
                         value_col   = "JI",
                         legend_name = "Jaccard\nIndex",
                         show_legend = TRUE) {
  
  ji_long$TF_x  <- factor(ji_long$TF_x, levels = tf_vec)
  ji_long$TF_y  <- factor(ji_long$TF_y, levels = tf_vec)
  ji_long$fill_ <- ji_long[[value_col]]
  
  ggplot(ji_long, aes(x = TF_x, y = TF_y, fill = fill_)) +
    geom_tile(color = "white", linewidth = 0.3) +
    # Tile labels; flip to white text when tile is dark
    geom_text(aes(label  = ifelse(is.na(fill_), "", sprintf("%.2f", fill_)),
                  colour = fill_ > 0.6),
              size = 2.8, show.legend = FALSE) +
    scale_colour_manual(values = c("TRUE" = "white", "FALSE" = "black")) +
    scale_fill_gradientn(
      # colours = c("#FCFDBF", "#FE9F6D", "#DE4968", "#8C2981", "#3B0F70"),
      colours = c("#F9F6F4", 
                  # "#F1DEDA", 
                  "#C684AA", "#6B4E80", "#4D4665"),
      limits  = c(0, 1),
      name    = legend_name
    ) +
    coord_fixed() +
    labs(title = title, x = "Transcription factor", y = NULL) +
    theme_classic(base_size = 10) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 10),
      axis.text.x     = element_text(angle = 45, hjust = 1),
      axis.ticks      = element_blank(),
      panel.border    = element_blank(),
      legend.position = if (show_legend) "right" else "none"
    )
}


# ── Helper: assemble interleaved panel (rows = CT, cols = gene | region) ──────
build_panel <- function(gplots, rplots, ct_vec) {
  
  interleaved <- unlist(
    lapply(ct_vec, function(ct) list(gplots[[ct]], rplots[[ct]])),
    recursive = FALSE
  )
  
  for (i in seq_along(ct_vec)) {
    interleaved[[2*i - 1]] <- interleaved[[2*i - 1]] +
      ggtitle(paste0(ct_vec[i], "\nTarget genes"))
    interleaved[[2*i    ]] <- interleaved[[2*i    ]] +
      ggtitle(paste0(ct_vec[i], "\nTarget regions"))
  }
  
  Reduce(`+`, interleaved) +
    plot_layout(ncol = 2, guides = "collect") &
    theme(legend.position = "right")
}


# ── Main ──────────────────────────────────────────────────────────────────────

# Read all requested sheets
message("Reading Excel sheets ...")
sheets_raw <- lapply(CT.v, function(ct) read_xlsx(xlsx_path, sheet = ct))
names(sheets_raw) <- CT.v

# Compute raw JI for all CTs
message("Computing pairwise JI ...")
all_ji_genes   <- lapply(CT.v, function(ct)
  compute_JI(sheets_raw[[ct]], TF.v, col = "Gene")   %>% mutate(CT = ct))
all_ji_regions <- lapply(CT.v, function(ct)
  compute_JI(sheets_raw[[ct]], TF.v, col = "Region") %>% mutate(CT = ct))

# For each TF pair, find max JI across all CTs → normalisation denominator
max_JI_genes <- bind_rows(all_ji_genes) %>%
  group_by(TF_x, TF_y) %>%
  summarise(JI_max = max(JI, na.rm = TRUE), .groups = "drop")

max_JI_regions <- bind_rows(all_ji_regions) %>%
  group_by(TF_x, TF_y) %>%
  summarise(JI_max = max(JI, na.rm = TRUE), .groups = "drop")

# ── Raw JI plots ──────────────────────────────────────────────────────────────
gene_plots <- lapply(CT.v, function(ct) {
  ji_long <- all_ji_genes[[which(CT.v == ct)]]
  make_heatmap(ji_long, TF.v, title = ct, legend_name = "JI\ngene")
})
names(gene_plots) <- CT.v

region_plots <- lapply(CT.v, function(ct) {
  ji_long <- all_ji_regions[[which(CT.v == ct)]]
  make_heatmap(ji_long, TF.v, title = ct, legend_name = "JI\nregion")
})
names(region_plots) <- CT.v

# ── CT-normalised JI plots ────────────────────────────────────────────────────
# JI_norm = JI / max(JI across all CT.v) for each TF pair → values in [0, 1]
gene_plots_norm <- lapply(CT.v, function(ct) {
  ji_long <- all_ji_genes[[which(CT.v == ct)]] %>%
    left_join(max_JI_genes, by = c("TF_x", "TF_y")) %>%
    mutate(JI_norm = ifelse(JI_max == 0 | is.na(JI_max), NA_real_, JI / JI_max))
  make_heatmap(ji_long, TF.v, title = ct,
               value_col = "JI_norm", legend_name = "Normalised\nJI\ngene")
})
names(gene_plots_norm) <- CT.v

region_plots_norm <- lapply(CT.v, function(ct) {
  ji_long <- all_ji_regions[[which(CT.v == ct)]] %>%
    left_join(max_JI_regions, by = c("TF_x", "TF_y")) %>%
    mutate(JI_norm = ifelse(JI_max == 0 | is.na(JI_max), NA_real_, JI / JI_max))
  make_heatmap(ji_long, TF.v, title = ct,
               value_col = "JI_norm", legend_name = "Normalised\nJI\nregion")
})
names(region_plots_norm) <- CT.v

# ── Assemble and save ─────────────────────────────────────────────────────────
panel_raw  <- build_panel(gene_plots,      region_plots,      CT.v)
panel_norm <- build_panel(gene_plots_norm, region_plots_norm, CT.v)

message("Saving figures ...")
ggsave(out_pdf, 
       panel_raw,  width = fig_width, height = fig_height,
       device = "pdf", useDingbats = FALSE)
ggsave(out_pdf_norm, 
       panel_norm, width = fig_width, height = fig_height,
       device = "pdf", useDingbats = FALSE)

message("Done!  Output written to:\n  ", out_pdf, "\n  ", out_pdf_norm)

