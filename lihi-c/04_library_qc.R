#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# Read the HiCUP output root and the directory for QC tables and figures.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: 04_library_qc.R <hicup_runs_dir> <output_dir> [16216 17759]")
}
base_dir <- normalizePath(args[1], mustWork = TRUE)
outdir <- args[2]
samples <- if (length(args) > 2L) args[3:length(args)] else {
  c("16216", "17759")
}
if (length(samples) != 2L || !setequal(samples, c("16216", "17759"))) stop("Use both disomy libraries only")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Accept known HiCUP column spellings and require nonnegative numeric counts.
sum_column <- function(dt, candidates, file) {
  hit <- candidates[candidates %in% names(dt)][1]
  if (is.na(hit)) stop("Missing column ", paste(candidates, collapse = "/"), " in ", file)
  values <- as.numeric(dt[[hit]])
  if (any(!is.finite(values)) || any(values < 0)) stop("Invalid counts in ", file)
  sum(values)
}

read_many <- function(files) {
  if (!length(files)) return(NULL)
  rbindlist(lapply(files, fread), fill = TRUE, use.names = TRUE)
}

# Read each donor's original filtering and duplicate-removal summary files.
read_hicup_qc <- function(sample) {
  sample_dir <- file.path(base_dir, sample, "output")
  filter_files <- Sys.glob(file.path(sample_dir, "hicup_filter_summary_*.txt"))
  dedup_files <- Sys.glob(file.path(sample_dir, "hicup_deduplicator_summary_*.txt"))
  if (length(filter_files) != 1L || length(dedup_files) != 1L) {
    stop("Expected exactly one filter and one deduplicator summary for ", sample, " in ", sample_dir)
  }
  filter_dt <- read_many(filter_files)
  dedup_dt <- read_many(dedup_files)

  # Total_pairs counts pairs entering HiCUP filtering, not all raw sequenced pairs.
  # Valid_pairs is before duplicate removal; unique_pairs is after duplicate removal.
  total_pairs <- sum_column(filter_dt, c("Total_pairs"), filter_files[1])
  valid_pairs <- sum_column(filter_dt, c("Valid_pairs"), filter_files[1])
  unique_pairs <- sum_column(dedup_dt, c("Unique_di-tags", "Unique_ditags"), dedup_files[1])
  cis_short <- sum_column(dedup_dt, c("Cis_<10kbp_of_uniques"), dedup_files[1])
  cis_long <- sum_column(dedup_dt, c("Cis_>10kbp_of_uniques"), dedup_files[1])
  trans <- sum_column(dedup_dt, c("Trans_of_uniques"), dedup_files[1])

  if (total_pairs <= 0 || valid_pairs <= 0 || unique_pairs <= 0 ||
      unique_pairs > valid_pairs || valid_pairs > total_pairs ||
      cis_short + cis_long + trans != unique_pairs) stop("Inconsistent HiCUP counts for ",sample)
  # Use the appropriate denominator for each fraction. Cis/trans fractions use uniques.
  data.table(
    Sample = sample,
    Total_pairs = total_pairs,
    Valid_pairs = valid_pairs,
    Unique_pairs = unique_pairs,
    Valid_pair_fraction = valid_pairs / total_pairs,
    Unique_pair_fraction = unique_pairs / valid_pairs,
    Cis_short_fraction = cis_short / unique_pairs,
    Cis_long_fraction = cis_long / unique_pairs,
    Trans_fraction = trans / unique_pairs,
    Duplicate_fraction = 1 - unique_pairs / valid_pairs,
    Cis_trans_ratio = (cis_short + cis_long) / trans
  )
}

# Save the numerical results before drawing the comparison plot.
qc <- rbindlist(lapply(sort(samples), read_hicup_qc))
fwrite(qc, file.path(outdir, "HiC_QC_metrics.tsv"), sep = "\t")

# Reshape the metrics so both libraries can be shown side by side.
plot_dt <- melt(qc, id.vars = "Sample", measure.vars = c(
  "Valid_pair_fraction", "Unique_pair_fraction", "Cis_short_fraction",
  "Cis_long_fraction", "Trans_fraction"
), variable.name = "Metric", value.name = "Fraction")
plot_dt[, Metric := factor(Metric,
  levels = c("Valid_pair_fraction", "Unique_pair_fraction", "Cis_short_fraction",
             "Cis_long_fraction", "Trans_fraction"),
  labels = c("Valid pairs\n/ total pairs", "Unique pairs\n/ valid pairs",
             "Cis <10 kb\n/ unique pairs", "Cis >10 kb\n/ unique pairs",
             "Trans\n/ unique pairs"))]

p <- ggplot(plot_dt, aes(Metric, 100 * Fraction, fill = Sample)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * Fraction)),
            position = position_dodge(width = 0.8), vjust = -0.35, size = 3) +
  scale_y_continuous("Read pairs (%)", limits = c(0, 105), breaks = seq(0, 100, 20)) +
  labs(x = NULL, fill = "Library", subtitle = "Hi-C library quality-control metrics") +
  theme_classic(base_size = 11) +
  theme(legend.position = "top")
ggsave(file.path(outdir, "HiC_QC_panel.pdf"), p, width = 9, height = 4.8)
ggsave(file.path(outdir, "HiC_QC_panel.png"), p, width = 9, height = 4.8, dpi = 300)
message("[DONE] QC outputs written to: ", outdir)
