###############################
### Metaplot B Lineage TFs  ###
###############################

###-------------------Note-------------------###
###-------      Run with R-4.4.2      -------###
###------------------------------------------###

##################
# Initialization #
##################

setwd("~/")
source("./00.Initialization.R")
dir.create(paste0(FIG_DIR, "02.metaplot_BLineage_eRegulon"))
dir.create(paste0(RES_DIR, "02.metaplot_BLineage_eRegulon"))

# Load addtional packages
library(caret)
library(pROC)
library(future)
library(furrr)
library(readxl)
library(soGGi)

##################
# Set parameters #
##################

xlsx_path <- "./data/CT_eGRN_DAR_Prune_Metadata.xlsx"   # per-cell-type eRegulon metadata
raw_xlsx_path <- "./data/eRegulon_Raw_Metadata.xlsx"    # global raw eRegulon metadata

# Load raw eRegulons
message("Reading raw eRegulon (+/+ only) ...")
raw_df <- read_xlsx(raw_xlsx_path, sheet = "+/+")

# Get cell type
CT.v <- c("HSC", "PreProB", "ProB-1", "ProB-2", "Large-PreB", "Small-PreB", "IM-B")
# Construct BAM paths following the existing naming convention
bam_dir <- "/work/DevM_analysis/data/bam_per_cell"

# Metaplots
TF.v <- c("EBF1", "PAX5", "FOXO1", "TCF3")

# Prepare TF binding regions from eRegulons
tf.idx <- TF.v[4]
TF_eRegulon.gr <- raw_df %>% 
  filter(TF == tf.idx) %>%
  mutate(
    chr = str_extract(Region, "^chr[^:]+"),
    start = as.integer(str_extract(Region, "(?<=:)[0-9]+")),
    end = as.integer(str_extract(Region, "(?<=-)[0-9]+"))
  ) %>%
  select(chr, start, end) %>%
  as(., "GRanges") %>%
  sort()

# Save to BED files
library(rtracklayer)
export(TF_eRegulon.gr, paste0("./results/02.metaplot_BLineage_eRegulon/", tf.idx, ".bed"), format = "bed")

#####################
# Generate metaplot #
#####################

# Compute profile for each CT separately
TF_profiles.list <- lapply(CT.v, function(ct) {
  bam <- file.path(bam_dir, ct, paste0(ct, ".dedup.bam"))
  message("Processing: ", ct)
  regionPlot(bam, TF_eRegulon.gr, style = "point", format = "bam")
})
names(TF_profiles.list) <- CT.v

signal_list <- lapply(CT.v, function(ct) {
  mat <- assays(TF_profiles.list[[ct]])[[1]]
  
  pos <- as.numeric(gsub("Point_Centre", "", colnames(mat)))
  
  total_reads = TF_profiles.list[[ct]]@metadata$AlignedReadsInBam
  
  data.frame(
    position = pos,
    coverage = colMeans(mat, na.rm = TRUE) / total_reads * 1e6,
    CT = ct
  )
})

# Organize plot df
signal_df <- bind_rows(signal_list)
## CT order
signal_df$CT <- factor(signal_df$CT, levels = intersect(CT_ORDER, unique(signal_df$CT)))
## Label data: take the rightmost position point for each CT
label_df <- signal_df %>%
  group_by(CT) %>%
  filter(position >= -200 & position <= 200) %>%   # restrict to peak region
  slice_max(coverage, n = 1) %>%
  ungroup()

# Save combined metaplot
Cairo::CairoPDF(paste0(FIG_DIR, "02.metaplot_BLineage_eRegulon/", tf.idx, "_B_Lineage.pdf"),
                width = 8, height = 6, family = "Arial")
ggplot(signal_df, aes(x = position, y = coverage, colour = CT)) +
  geom_line(linewidth = 1, na.rm = TRUE) +
  geom_text_repel(
    data           = label_df,
    aes(label      = CT),
    size           = 3,
    fontface       = "bold",
    direction      = "y",
    nudge_x        = 500,
    nudge_y        = 0.5,
    segment.size   = 0.3,
    segment.colour = "grey50",
    show.legend    = FALSE,
    max.overlaps = Inf
  ) +
  scale_color_manual(values = COL) +
  labs(x = "Distance from center (bp)", y = "Mean coverage (RPM)",
       title = paste0("Coverage Plot for Regions from ", tf.idx, " eRegulon"), 
       colour = NULL) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = .5, colour = "black", size = 17))
dev.off()

