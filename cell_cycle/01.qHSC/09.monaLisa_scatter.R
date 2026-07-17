# ============================================================================
#  monaLisa binned motif enrichment: Up vs Dn, stratified by genomic region
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
  library(JASPAR2024)
})

# ---- 0. USER PARAMETERS ----------------------------------------------------
setwd("~/")
res_dir <- "./results/09.monaLisa_scatter/"
fig_dir <- "./plots/09.monaLisa_scatter/"
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

## Load Lisha qHSC/cHSC regions
diff_peaks.df <- read.table(
  "~/21.G0_vs_Cyc_diff_peak.bed",
  header = TRUE, sep = "\t"
)
### Filter by FDR threshold
diff_peaks.df <- diff_peaks.df %>%
  filter(FDR <= 0.05)
### Split into Up (more accessible in G0/qHSC) and Dn
RegionSets_Up.df <- diff_peaks.df %>% filter(Log2FC > 0)
RegionSets_dn.df <- diff_peaks.df %>% filter(Log2FC < 0)
### Info
message("Up-in-G0: ", nrow(RegionSets_Up.df),
        "  Dn-in-G0: ", nrow(RegionSets_dn.df))

## Build metadata
RegionSets_Up.metadata <- distinct(RegionSets_Up.df, seqnames, start, end)
RegionSets_dn.metadata <- distinct(RegionSets_dn.df, seqnames, start, end)
### Info
message("Up metadata: ", nrow(RegionSets_Up.metadata),
        "  Dn metadata: ", nrow(RegionSets_dn.metadata))

## Genome setup
txdb   <- TxDb.Hsapiens.UCSC.hg38.knownGene
genome <- BSgenome.Hsapiens.UCSC.hg38

# ---- 1. READ PEAKS ---------------------------------------------------------
up_peaks <- makeGRangesFromDataFrame(RegionSets_Up.metadata) |> sort()
dn_peaks <- makeGRangesFromDataFrame(RegionSets_dn.metadata) |> sort()
message("Up peaks: ", length(up_peaks), "   Dn peaks: ", length(dn_peaks))

# ---- 2. CHIPSEEKER ANNOTATION ---------------------------------------------
## Annotate peaks
anno_up <- annotatePeak(up_peaks, tssRegion = c(-3000, 3000),
                        TxDb = txdb, level = "transcript",
                        verbose = FALSE)
anno_dn <- annotatePeak(dn_peaks, tssRegion = c(-3000, 3000),
                        TxDb = txdb, level = "transcript",
                        verbose = FALSE)

gr_up <- as.GRanges(anno_up)
gr_dn <- as.GRanges(anno_dn)
### Save Annotation
saveRDS(list(anno_up = anno_up,
             anno_dn = anno_dn,
             gr_up = gr_up,
             gr_dn = gr_dn), file = paste0(res_dir, "up_dn_peak_annotation.rds"))

## Subset promoters & distal regions
### ChIPseeker writes things like "Promoter (<=1kb)", "Distal Intergenic". Match by prefix.
region_mask <- function(gr, region_label) startsWith(gr$annotation, region_label)
subset_region <- function(region_label) {
  up_sel <- gr_up[region_mask(gr_up, region_label)]
  dn_sel <- gr_dn[region_mask(gr_dn, region_label)]
  gr <- c(up_sel, dn_sel)
  bins <- factor(c(rep("dn", length(dn_sel)),
                   rep("up", length(up_sel))),
                 levels = c("dn", "up"))
  list(gr = gr, bins = bins)
}
promoter_data <- subset_region("Promoter (<=1kb)")
distal_data <- subset_region("Distal Intergenic")
### No-subset all peaks
all_data <- list(
  gr   = c(gr_dn, gr_up),
  bins = factor(c(rep("dn", length(gr_dn)), rep("up", length(gr_up))),
                levels = c("dn", "up"))
)
### Info
message(sprintf("Promoter (<=1kb):   dn=%f  up=%f",
                sum(promoter_data$bins == "dn")/length(dn_peaks),
                sum(promoter_data$bins == "up")/length(up_peaks)))
message(sprintf("Distal Intergenic: dn=%f  up=%f",
                sum(distal_data$bins == "dn")/length(dn_peaks),
                sum(distal_data$bins == "up")/length(up_peaks)))

# ---- 3. MOTIFS (JASPAR CORE vertebrates) ----------------------------------
jdb  <- JASPAR2024::JASPAR2024()
con  <- RSQLite::dbConnect(RSQLite::SQLite(), jdb@db)
pwms <- TFBSTools::getMatrixSet(con,
                                opts = list(tax_group = "vertebrates", collection = "CORE", matrixtype = "PWM"))
RSQLite::dbDisconnect(con)
message("Motifs loaded: ", length(pwms))

# ---- 4. PREPARE SEQUENCES & RUN monaLisa ----------------------------------
peak_width    <- 500
prep_and_run <- function(region_data, tag, res_dir) {
  gr   <- region_data$gr
  bins <- region_data$bins
  
  # Resize to fixed width centered on peak center
  gr  <- trim(suppressWarnings(resize(gr, width = peak_width, fix = "center")))
  keep <- width(gr) == peak_width &
    as.character(seqnames(gr)) %in% seqlevels(genome)
  gr   <- gr[keep]
  bins <- droplevels(bins[keep])
  
  # Drop seqlevels that are gone after filtering
  seqlevels(gr) <- seqlevelsInUse(gr)
  
  # Extract sequences
  seqs <- BSgenome::getSeq(genome, gr)
  
  message(sprintf("[%s] scanning %d peaks across %d motifs...",
                  tag, length(gr), length(pwms)))
  
  se <- calcBinnedMotifEnrR(
    seqs       = seqs,
    bins       = bins,
    pwmL       = pwms,
    background = "otherBins",
    # BPPARAM    = BiocParallel::SerialParam(),
    BPPARAM = MulticoreParam(workers = 10),
    verbose    = TRUE
  )
  
  saveRDS(se, file.path(res_dir, sprintf("monalisa_%s.rds", tag)))
  se
}

se_promoter <- prep_and_run(promoter_data, "Promoter_leq1kb", res_dir)
se_distal <- prep_and_run(distal_data,   "Distal_Intergenic", res_dir)
se_all <- prep_and_run(all_data, "All", res_dir)

# ---- 5. HEATMAP supplement to Fig.3e -------------------------------------------
plot_fig3e_style <- function(se, tag, fig_dir, sig_cutoff, enr_cutoff) {
  # Filter: significant + |enr| above threshold in at least one bin
  sel <- apply(assay(se, "negLog10Padj"), 1,
               function(x) any(x > sig_cutoff, na.rm = TRUE)) &
    apply(assay(se, "log2enr"), 1,
          function(x) any(abs(x) > enr_cutoff, na.rm = TRUE))
  se_sel <- se[sel, ]
  message(sprintf("[%s] %d motifs pass filter", tag, nrow(se_sel)))
  if (nrow(se_sel) == 0) return(invisible(NULL))
  
  pdf_path <- paste0(fig_dir, sprintf("monalisa_heatmap_%s.pdf", tag))
  pdf(pdf_path, width = 10, height = max(4, nrow(se_sel) * 0.2))
  plotMotifHeatmaps(
    x              = se_sel,
    which.plots    = c("log2enr", "negLog10Padj"),
    width          = 2.5,
    width.seqlogo  = 1.2,
    cluster        = TRUE,
    maxEnr         = 2,
    maxSig         = 10,
    show_motif_GC  = TRUE,
    show_dendrogram = FALSE,
    show_seqlogo   = TRUE
  )
  dev.off()
  message("  wrote ", pdf_path)
}

plot_fig3e_style(se_promoter, "Promoter_leq1kb", fig_dir, sig_cutoff = 3,enr_cutoff)
plot_fig3e_style(se_distal,   "Distal_Intergenic", fig_dir, sig_cutoff = 3, enr_cutoff)
plot_fig3e_style(se_all,   "All", fig_dir, sig_cutoff = 4, enr_cutoff)

# ---- 6. Scatter: Promoter vs Distal -----------------------------------------
sig_cutoff     <- 3   # for the prom/dist filter (which TFs are plotted)
sig_cutoff_all <- 4   # for selecting TFs to label (TFs that also pass se_all)
filter_mode    <- "any"   # "any" | "both"

if (!exists("se_promoter")) se_promoter <- readRDS(file.path(res_dir, "monalisa_Promoter_leq1kb.rds"))
if (!exists("se_distal"))   se_distal   <- readRDS(file.path(res_dir, "monalisa_Distal_Intergenic.rds"))
if (!exists("se_all"))      se_all      <- readRDS(file.path(res_dir, "monalisa_All.rds"))

# log2enr[,"up"] is the Up-vs-Dn enrichment; for -log10 adj.p take max across bins
tidy_se <- function(se) {
  l2e  <- assay(se, "log2enr")
  padj <- assay(se, "negLog10Padj")
  nm   <- rowData(se)$motif.name %||% rownames(se)
  data.frame(motif_id   = rownames(se),
             motif_name = as.character(nm),
             log2enr    = l2e[, "up"],
             nlp        = pmax(padj[, "up"], padj[, "dn"], na.rm = TRUE))
}

`%||%` <- function(a, b) if (is.null(a)) b else a   # tiny helper
prom <- tidy_se(se_promoter); colnames(prom)[3:4] <- c("log2enr_prom", "nlp_prom")
dist <- tidy_se(se_distal);   colnames(dist)[3:4] <- c("log2enr_dist", "nlp_dist")

tbl <- merge(prom, dist, by = c("motif_id", "motif_name")) %>%
  filter(!is.na(log2enr_prom), !is.na(log2enr_dist))
tbl$nlp_mean <- (tbl$nlp_prom + tbl$nlp_dist) / 2

# significance per region; direction from sign(log2enr)
sig_prom <- tbl$nlp_prom >= sig_cutoff
sig_dist <- tbl$nlp_dist >= sig_cutoff
up_prom  <- tbl$log2enr_prom > 0
up_dist  <- tbl$log2enr_dist > 0

tbl$class <- dplyr::case_when(
  sig_prom &  sig_dist &  up_prom &  up_dist        ~ "Enriched in both for qHSC",
  sig_prom &  sig_dist & !up_prom & !up_dist        ~ "Enriched in both for non-qHSC",
  sig_prom &  sig_dist & (up_prom != up_dist)       ~ "Opposite",
  sig_prom & !sig_dist                              ~ "Promoter-specific",
  !sig_prom &  sig_dist                              ~ "Distal-specific",
  TRUE                                               ~ "Other"
)
tbl$kept <- if (filter_mode == "both") sig_prom & sig_dist else sig_prom | sig_dist

# TFs that also pass se_all -> these get labels
ids_in_all <- rownames(se_all)[
  apply(assay(se_all, "negLog10Padj"), 1, function(x) any(x > sig_cutoff_all, na.rm = TRUE))
]
tbl$in_all <- tbl$motif_id %in% ids_in_all

write.csv(tbl, file.path(res_dir, "TF_promoter_vs_distal_table.csv"), row.names = FALSE)

# keep only classes actually present, in a fixed legend order
class_levels <- c("Enriched in both for qHSC",
                  "Enriched in both for non-qHSC",
                  "Promoter-specific",
                  "Distal-specific",
                  "Opposite")
class_colors <- c("Enriched in both for qHSC"     = "#E41A1C",
                  "Enriched in both for non-qHSC" = "#377EB8",
                  "Promoter-specific"             = "#FF7F00",
                  "Distal-specific"               = "#E41A1C",
                  "Opposite"                      = "#4DAF4A")

tbl_plot <- tbl[tbl$kept, ]
tbl_plot$class <- factor(tbl_plot$class, levels = class_levels)
tbl_plot <- tbl_plot[order(tbl_plot$in_all), ]   # labeled points on top

message("Class breakdown:"); print(table(tbl_plot$class))

lim <- max(abs(c(tbl_plot$log2enr_prom, tbl_plot$log2enr_dist)), na.rm = TRUE) * 1.1

p <- ggplot(tbl_plot, aes(log2enr_prom, log2enr_dist)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey40") +
  geom_point(aes(size = nlp_mean, fill = class),
             shape = 21, color = "black", stroke = 0.3, alpha = 0.9) +
  ggrepel::geom_text_repel(
    data = subset(tbl_plot, in_all),
    aes(label = motif_name, color = class),
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
       subtitle = sprintf("%d motifs (-log10 adj.p > %g)",
                          nrow(tbl_plot), sig_cutoff)) +
  guides(fill = guide_legend(override.aes = list(size = 5),
                             title.position = "top", ncol = 1),
         size = guide_legend(title.position = "top", ncol = 1)) +
  theme_classic(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold", hjust = 0, size = 13),
        plot.subtitle    = element_text(face = "bold", hjust = 0, size = 9),
        legend.position  = "bottom",
        legend.box       = "horizontal")

ggsave(file.path(fig_dir, "scatter_promoter_vs_distal.pdf"),
       p, width = 6, height = 7.5, device = cairo_pdf)
