############################################################
# Hi-C replicate QC plots from HiCDCPlus results
# Observed-count analyses only
############################################################

library(data.table)
library(ggplot2)
library(scales)

# ----------------------------------------------------------
# Parameters
# ----------------------------------------------------------

# Supply donor 16216 first and donor 17759 second, using their 25-kb call tables.
# This script writes the numerical QC tables and the corresponding PDF plots.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<3L) stop("Usage: Rscript 05_replicate_qc.R rep1_25kb.txt.gz rep2_25kb.txt.gz outdir [rep1_name=16216] [rep2_name=17759]")
f_rep1 <- normalizePath(args[1],mustWork=TRUE)
f_rep2 <- normalizePath(args[2],mustWork=TRUE)
outdir <- args[3]
rep1_name <- if(length(args)>=4L) args[4] else "16216"
rep2_name <- if(length(args)>=5L) args[5] else "17759"
if(!identical(c(rep1_name,rep2_name),c("16216","17759"))) stop("Supply 16216 first and 17759 second")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

# Include cis interactions from one 25-kb bin separation through 5 Mb.
min_distance <- 25000
max_distance <- 5000000

# At most this many rows are used to draw the density plot.
# Correlations are calculated using all eligible interactions.
plot_n <- 500000

# Minimum number of interactions needed to calculate a correlation
# within a genomic-distance stratum.
min_interactions_per_bin <- 100

# Distance strata used for Panel B.
distance_breaks <- c(
  25000,
  50000,
  100000,
  200000,
  500000,
  1000000,
  2000000,
  5000000 + 1
)

distance_labels <- c(
  "25-50 kb",
  "50-100 kb",
  "100-200 kb",
  "200-500 kb",
  "500 kb-1 Mb",
  "1-2 Mb",
  "2-5 Mb"
)

set.seed(1)

# ----------------------------------------------------------
# Load only the columns needed
# ----------------------------------------------------------

message("[INFO] Loading replicate 1")

rep1 <- fread(
  f_rep1,
  select = c(
    "chrI", "startI", "endI",
    "chrJ", "startJ", "endJ",
    "D", "counts", "pvalue", "qvalue"
  ),
  showProgress = TRUE
)

message("[INFO] Loading replicate 2")

rep2 <- fread(
  f_rep2,
  select = c(
    "chrI", "startI", "endI",
    "chrJ", "startJ", "endJ",
    "D", "counts", "pvalue", "qvalue"
  ),
  showProgress = TRUE
)

message("[INFO] Replicate 1 rows: ", format(nrow(rep1), big.mark = ","))
message("[INFO] Replicate 2 rows: ", format(nrow(rep2), big.mark = ","))

# Check that each input has unique bin pairs and the expected bin width.
for (z in list(rep1, rep2)) {
 if (anyDuplicated(z[,.(chrI,startI,endI,chrJ,startJ,endJ)])) stop("Duplicate coordinates in input table")
 if (!median(z$endI-z$startI) %in% c(24999,25000) || !median(z$endJ-z$startJ) %in% c(24999,25000)) stop("Use 25-kb interaction tables")
}
# ----------------------------------------------------------
# Combine replicate tables
# ----------------------------------------------------------

coord_cols <- c(
  "chrI", "startI", "endI",
  "chrJ", "startJ", "endJ",
  "D"
)

rows_match <- (
  nrow(rep1) == nrow(rep2) &&
    identical(rep1[, ..coord_cols], rep2[, ..coord_cols])
)

message("[INFO] Rows are in identical order: ", rows_match)

if (rows_match) {
  message("[INFO] Combining by row position")

  merged <- rep1[, .(
    chrI, startI, endI,
    chrJ, startJ, endJ,
    D,
    counts_rep1 = counts,
    pvalue_rep1 = pvalue,
    qvalue_rep1 = qvalue
  )]

  merged[, `:=`(
    counts_rep2 = rep2$counts,
    pvalue_rep2 = rep2$pvalue,
    qvalue_rep2 = rep2$qvalue
  )]

} else {
  message("[INFO] Rows differ; merging by coordinates")

  setnames(
    rep1,
    c("counts", "pvalue", "qvalue"),
    c("counts_rep1", "pvalue_rep1", "qvalue_rep1")
  )

  setnames(
    rep2,
    c("counts", "pvalue", "qvalue"),
    c("counts_rep2", "pvalue_rep2", "qvalue_rep2")
  )

  merged <- merge(
    rep1,
    rep2,
    by = coord_cols,
    all = FALSE,
    sort = FALSE
  )
}

rm(rep1, rep2)
gc()

message("[INFO] Combined rows: ", format(nrow(merged), big.mark = ","))

# ----------------------------------------------------------
# Eligible interactions for reproducibility analysis
# ----------------------------------------------------------

# Exclude D = 0 and exclude rows where both replicates have zero counts.
# Shared zeroes can otherwise inflate apparent replicate concordance.
count_qc <- merged[
  chrI == chrJ &
    D >= min_distance &
    D <= max_distance &
    (counts_rep1 > 0 | counts_rep2 > 0)
]

message(
  "[INFO] Eligible nonzero-in-at-least-one-replicate interactions between ",
  comma(min_distance),
  " and ",
  comma(max_distance),
  " bp: ",
  comma(nrow(count_qc))
)

if (nrow(count_qc) == 0) {
  stop("No eligible interactions remained after filtering.")
}

############################################################
# PANEL B: Observed-count concordance across genomic distance
############################################################

count_qc[, Distance_bin := cut(
  D,
  breaks = distance_breaks,
  labels = distance_labels,
  include.lowest = TRUE,
  right = FALSE
)]

count_qc <- count_qc[!is.na(Distance_bin)]

safe_cor <- function(x, y, method) {
  keep <- complete.cases(x, y)

  if (sum(keep) < min_interactions_per_bin) {
    return(NA_real_)
  }

  # cor() is undefined when either vector has zero variance.
  if (sd(x[keep]) == 0 || sd(y[keep]) == 0) {
    return(NA_real_)
  }

  cor(x[keep], y[keep], method = method)
}

distance_correlations <- count_qc[, .(
  N = .N,
  Pearson = safe_cor(
    log1p(counts_rep1),
    log1p(counts_rep2),
    method = "pearson"
  ),
  Spearman = safe_cor(
    counts_rep1,
    counts_rep2,
    method = "spearman"
  )
), by = Distance_bin]

distance_correlations[, Distance_bin := factor(
  Distance_bin,
  levels = distance_labels
)]

setorder(distance_correlations, Distance_bin)

message("\n[INFO] Distance-stratified correlations:")
print(distance_correlations)

fwrite(
  distance_correlations,
  file.path(outdir, "panel_B_observed_count_correlations_by_distance.tsv"),
  sep = "\t"
)

# Plot Pearson as the primary metric; retain Spearman in the TSV.
p_distance_counts <- ggplot(
  distance_correlations,
  aes(
    x = Distance_bin,
    y = Pearson,
    group = 1
  )
) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_point(size = 2.5, na.rm = TRUE) +
  geom_text(
    aes(label = paste0("n=", comma(N))),
    vjust = -0.8,
    size = 3,
    na.rm = TRUE
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0.02, 0.12))
  ) +
  labs(
    subtitle = "Replicate concordance across genomic distance",
    x = "Genomic distance",
    y = "Pearson correlation of log1p contact counts"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(face = "bold"),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

ggsave(
  file.path(outdir, "panel_B_count_correlation_by_distance.pdf"),
  p_distance_counts,
  width = 6,
  height = 4.5
)

############################################################
# PANEL C: Overall observed-count replicate concordance
############################################################

pearson_counts <- cor(
  log1p(count_qc$counts_rep1),
  log1p(count_qc$counts_rep2),
  method = "pearson",
  use = "complete.obs"
)

spearman_counts <- cor(
  count_qc$counts_rep1,
  count_qc$counts_rep2,
  method = "spearman",
  use = "complete.obs"
)

message("[INFO] Overall count Pearson correlation: ", round(pearson_counts, 3))
message("[INFO] Overall count Spearman correlation: ", round(spearman_counts, 3))

# Sample rows only for drawing. Correlations above use every eligible row.
if (nrow(count_qc) > plot_n) {
  count_plot_dt <- count_qc[sample.int(.N, plot_n)]
} else {
  count_plot_dt <- copy(count_qc)
}

count_annotation <- sprintf(
  "Pearson r = %.3f\nSpearman rho = %.3f\nn = %s",
  pearson_counts,
  spearman_counts,
  comma(nrow(count_qc))
)

p_counts <- ggplot(
  count_plot_dt,
  aes(
    x = log10(counts_rep1 + 1),
    y = log10(counts_rep2 + 1)
  )
) +
  geom_bin_2d(bins = 80) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = count_annotation,
    hjust = 1.05,
    vjust = -0.4,
    size = 3.5
  ) +
  labs(
    subtitle = "Observed contact-count concordance",
    x = paste0(rep1_name, ": log10(count + 1)"),
    y = paste0(rep2_name, ": log10(count + 1)"),
    fill = "Number of\ninteractions"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(
  file.path(outdir, "panel_C_count_concordance.pdf"),
  p_counts,
  width = 5.5,
  height = 4.5
)

############################################################
# Additional summary: significant-interaction overlap
############################################################

# Calculate significance overlap over the same distance range, but do not
# require a nonzero observed count because significance itself defines this
# secondary comparison.
sig_qc <- merged[
  chrI == chrJ &
    D >= min_distance &
    D <= max_distance
]

sig_qc[, significant_rep1 := !is.na(qvalue_rep1) & qvalue_rep1 <= 0.05]
sig_qc[, significant_rep2 := !is.na(qvalue_rep2) & qvalue_rep2 <= 0.05]

significance_summary <- sig_qc[, .(
  total_tested = .N,
  significant_rep1 = sum(significant_rep1),
  significant_rep2 = sum(significant_rep2),
  significant_both = sum(significant_rep1 & significant_rep2),
  significant_either = sum(significant_rep1 | significant_rep2)
)]

significance_summary[, overlap_rep1 := fifelse(
  significant_rep1 > 0,
  significant_both / significant_rep1,
  NA_real_
)]

significance_summary[, overlap_rep2 := fifelse(
  significant_rep2 > 0,
  significant_both / significant_rep2,
  NA_real_
)]

significance_summary[, jaccard := fifelse(
  significant_either > 0,
  significant_both / significant_either,
  NA_real_
)]

message("\n[INFO] Significant-interaction overlap:")
print(significance_summary)

fwrite(
  significance_summary,
  file.path(outdir, "significant_interaction_overlap.tsv"),
  sep = "\t"
)

############################################################
# PANEL D: Significant-interaction overlap across distance
############################################################

# Distance strata used for Panel D.
distance_breaks <- c(
  25000,
  200000,
  500000,
  1000000,
  2000000,
  5000000 + 1
)

distance_labels <- c(
  "25-200 kb",
  "200-500 kb",
  "500 kb-1 Mb",
  "1-2 Mb",
  "2-5 Mb"
)

# Use the same genomic-distance strata as Panel B.
sig_qc[, Distance_bin := cut(
  D,
  breaks = distance_breaks,
  labels = distance_labels,
  include.lowest = TRUE,
  right = FALSE
)]

sig_qc <- sig_qc[!is.na(Distance_bin)]

# For each distance bin, form the 2 x 2 table:
#
#                         Replicate 2
#                      Significant  Not significant
# Replicate 1 Significant     both       rep1_only
#             Not significant rep2_only   neither
#
# Fisher's exact-test OR > 1 means significant interactions in one
# replicate are enriched for significance in the other replicate.
get_fisher_stats <- function(both, rep1_only, rep2_only, neither) {
  tab <- matrix(
    c(both, rep1_only, rep2_only, neither),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      Replicate1 = c("Significant", "Not_significant"),
      Replicate2 = c("Significant", "Not_significant")
    )
  )

  ft <- fisher.test(tab)

  list(
    fisher_OR = unname(ft$estimate),
    fisher_CI_low = unname(ft$conf.int[1]),
    fisher_CI_high = unname(ft$conf.int[2]),
    fisher_pvalue = ft$p.value
  )
}

significance_by_distance <- sig_qc[, {
  both <- sum(significant_rep1 & significant_rep2)
  rep1_only <- sum(significant_rep1 & !significant_rep2)
  rep2_only <- sum(!significant_rep1 & significant_rep2)
  neither <- sum(!significant_rep1 & !significant_rep2)

  fisher_stats <- get_fisher_stats(
    both = both,
    rep1_only = rep1_only,
    rep2_only = rep2_only,
    neither = neither
  )

  significant_rep1 <- both + rep1_only
  significant_rep2 <- both + rep2_only
  significant_either <- both + rep1_only + rep2_only

  c(
    list(
      total_tested = .N,
      significant_rep1 = significant_rep1,
      significant_rep2 = significant_rep2,
      significant_both = both,
      rep1_only = rep1_only,
      rep2_only = rep2_only,
      neither = neither,
      overlap_rep1 = if (significant_rep1 > 0) {
        both / significant_rep1
      } else {
        NA_real_
      },
      overlap_rep2 = if (significant_rep2 > 0) {
        both / significant_rep2
      } else {
        NA_real_
      },
      jaccard = if (significant_either > 0) {
        both / significant_either
      } else {
        NA_real_
      }
    ),
    fisher_stats
  )
}, by = Distance_bin]

significance_by_distance[, Distance_bin := factor(
  Distance_bin,
  levels = distance_labels
)]

setorder(significance_by_distance, Distance_bin)

# Benjamini-Hochberg correction across distance bins.
significance_by_distance[, fisher_qvalue := p.adjust(
  fisher_pvalue,
  method = "BH"
)]

message("\n[INFO] Significant-interaction overlap by distance:")
print(significance_by_distance)

fwrite(
  significance_by_distance,
  file.path(outdir, "panel_D_significant_overlap_by_distance.tsv"),
  sep = "\t"
)

# Plot Fisher OR on a log scale. The dashed horizontal line at OR = 1
# represents no enrichment of significance overlap.
# Points are labelled with the number significant in both replicates.
p_significance_OR <- ggplot(
  significance_by_distance,
  aes(
    x = Distance_bin,
    y = fisher_OR,
    group = 1
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed"
  ) +
  geom_errorbar(
    aes(
      ymin = fisher_CI_low,
      ymax = fisher_CI_high
    ),
    width = 0.15,
    na.rm = TRUE
  ) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_point(size = 2.5, na.rm = TRUE) +
  geom_text(
    aes(label = paste0("both=", comma(significant_both))),
    vjust = -0.8,
    size = 3,
    na.rm = TRUE
  ) +
  scale_y_log10(
    breaks = scales::log_breaks(n = 6),
    labels = label_number(accuracy = 0.1),
    expand = expansion(mult = c(0.08, 0.18))
  ) +
  labs(
    title = "D",
    subtitle = "Significant-interaction concordance across genomic distance",
    x = "Genomic distance",
    y = "Fisher's exact-test odds ratio (95% CI)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(face = "bold"),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

ggsave(
  file.path(outdir, "panel_D_significant_overlap_OR_by_distance.pdf"),
  p_significance_OR,
  width = 6,
  height = 4.5
)

# Optional companion plot showing the overlap fractions directly.
overlap_plot_dt <- melt(
  significance_by_distance,
  id.vars = c("Distance_bin", "significant_both"),
  measure.vars = c("overlap_rep1", "overlap_rep2", "jaccard"),
  variable.name = "Metric",
  value.name = "Fraction"
)

overlap_plot_dt[, Metric := factor(
  Metric,
  levels = c("overlap_rep1", "overlap_rep2", "jaccard"),
  labels = c(
    paste0("Shared / ", rep1_name, " significant"),
    paste0("Shared / ", rep2_name, " significant"),
    "Jaccard index"
  )
)]

p_significance_fraction <- ggplot(
  overlap_plot_dt,
  aes(
    x = Distance_bin,
    y = Fraction,
    group = Metric,
    linetype = Metric,
    shape = Metric
  )
) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_point(size = 2.5, na.rm = TRUE) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    labels = label_percent(accuracy = 1)
  ) +
  labs(
    subtitle = "Significant-interaction overlap across genomic distance",
    x = "Genomic distance",
    y = "Overlap",
    linetype = NULL,
    shape = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(face = "bold"),
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "top"
  )

ggsave(
  file.path(outdir, "panel_D_significant_overlap_fraction_by_distance.pdf"),
  p_significance_fraction,
  width = 6.5,
  height = 4.5
)

############################################################
# Save overall observed-count correlation statistics
############################################################

correlation_summary <- data.table(
  Comparison = "Observed counts",
  Pearson = pearson_counts,
  Spearman = spearman_counts,
  N = nrow(count_qc),
  Minimum_distance = min_distance,
  Maximum_distance = max_distance,
  Shared_zeroes_excluded = TRUE
)

fwrite(
  correlation_summary,
  file.path(outdir, "replicate_correlations.tsv"),
  sep = "\t"
)

message("\n[INFO] Wrote results to: ", outdir)

# Global overlap enrichment; distance-stratified tests and their BH q-values are above.
a <- sum(sig_qc$significant_rep1 & sig_qc$significant_rep2)
b <- sum(sig_qc$significant_rep1 & !sig_qc$significant_rep2)
c <- sum(!sig_qc$significant_rep1 & sig_qc$significant_rep2)
d <- sum(!sig_qc$significant_rep1 & !sig_qc$significant_rep2)
ft <- fisher.test(matrix(c(a,b,c,d),nrow=2))
fwrite(data.table(odds_ratio=unname(ft$estimate),ci_low=ft$conf.int[1],ci_high=ft$conf.int[2],pvalue=ft$p.value),file.path(outdir,"global_fisher.tsv"),sep="\t")
capture.output(sessionInfo(),file=file.path(outdir,"R_sessionInfo.txt"))
