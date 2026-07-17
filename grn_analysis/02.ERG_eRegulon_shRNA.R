#####################################
### ERG_eRegulon_shRNA validation ###
#####################################

###-------------------Note-------------------###
###-------      Run with R-4.4.2      -------###
###------------------------------------------###

##################
# Initialization #
##################

setwd("~/")
dir.create(paste0(FIG_DIR, "02.ERG_eRegulon_shRNA"))
dir.create(paste0(RES_DIR, "02.ERG_eRegulon_shRNA"))

################
# Set eRegulon #
################

raw_xlsx_path <- "./data/eRegulon_Raw_Metadata.xlsx"

# Load raw eRegulons
message("Reading raw eRegulon (+/+ only) ...")
raw_df <- readxl::read_xlsx(raw_xlsx_path, sheet = "+/+")
TG.v <- raw_df %>% 
  filter(TF == "ERG") %>%
  pull(Gene) %>%
  unique()

## ---- shERG barplot: ERG & NKAIN2 (style of panel D) ------------------------

library(readxl)
library(data.table)
library(ggplot2)

xlsx_path  <- "./data/GSE158795_summary_AS_CD34kd_deSeq_analysis.xlsx"
genes_plot <- c("ERG", "NKAIN2")
FLIP_TO_KD <- TRUE
thr        <- 0.2

raw <- as.data.table(read_excel(xlsx_path))
lfc_col <- grep("log2FoldChange", names(raw), value = TRUE)[1]
raw[, sym := sub("^\\d+_+", "", gene)]

bar <- raw[sym %in% genes_plot,
           .(gene = sym, l2fc = as.numeric(get(lfc_col)))]
if (FLIP_TO_KD) bar[, l2fc := -l2fc]
bar[, gene := factor(gene, levels = genes_plot)]

COL <- c(ERG = "#2E7D32", NKAIN2 = "#E07B39")

p_bar <- ggplot(bar, aes(gene, l2fc, fill = gene)) +
  geom_col(width = 0.7, colour = NA) +
  geom_hline(yintercept = 0, linewidth = .7, colour = "black") +
  scale_fill_manual(values = COL, guide = "none") +
  scale_y_continuous(limits = c(-0.7, 0.1), breaks = c(-0.7, 0),
                     expand = c(0, 0)) +
  labs(x = NULL, y = "Expression (log2 foldchange)",
       title = expression(shERG~"(CD34"^"+"*", RNAseq)")) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 12) +
  theme(plot.title  = element_text(hjust = 0.5, size = 12),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin  = margin(t = 6, r = 12, b = 10, l = 6),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1))

ggsave("./plots/02.ERG_eRegulon_shRNA/shERG_ERG_NKAIN2_barplot.pdf", p_bar, width = 2.6, height = 4.5,
       device = cairo_pdf)
print(bar)



