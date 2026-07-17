#------- Check AGM HSC Cell Cycle phase percentage in HM paper -------#

# Load env
setwd("~/")
source("./00.initialization.R")

# Load HM CS14&15 AGM data
Calvanese_EDFig1f_AGMHema.SeuratObj <- readRDS("../results/Calvanese_EDFig1f_Reprocess_SeuratObj.rds")
Calvanese_EDFig1f_AGMHema.SeuratObj$organ <- "AGM"

# Rename Celltype
cluster_to_CT <- c(
  "1" = "HSC",
  "0" = "Mo/Mφ",
  "3" = "Gr",
  "2" = "Mo/Mφ",
  "4" = "Mo/Mφ",
  "7" = "Ly",
  "6" = "Mo/Mφ",
  "5" = "Gr",
  "8" = "Mo/Mφ"
)
md <- Calvanese_EDFig1f_AGMHema.SeuratObj@meta.data
md$CT_HM <- cluster_to_CT[as.character(md$seurat_clusters)]
Calvanese_EDFig1f_AGMHema.SeuratObj@meta.data <- md
CT_COL <- c("HSC" = "#4DAF4A", "Mo/Mφ" = "#E41A1C", "Ly" = "#377EB8", "Gr" = "#4DD2FF")
CT_ORDER <- names(CT_COL)

# Check Cell type
table(Calvanese_EDFig1f_AGMHema.SeuratObj$CT_HM)
p <- DimPlot(Calvanese_EDFig1f_AGMHema.SeuratObj, group.by = "CT_HM",
             label = T, label.size = 7,
             cols = CT_COL) +
  ggtitle("UMAP with cell type \nin Calvanese et al. AGM CS14&15 haematopoietic cells") +
  theme(axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1),
        axis.text = element_text(face = "bold", colour = "black"),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = .5, colour = "black", size = 15))

ggsave(
  filename = "plots/HM_AGM_UMAP.pdf",
  p,
  width = 8, height = 7,
  device = cairo_pdf
)

# Get cell cycle genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
## Set color
CC_COL <- c(
  "G1"   = "#2a6da1",
  "S"    = "#E6762B",
  "G2M"  = "#329848"
)

# Find alias
alias_table.S <- resolve_gene_aliases(s.genes, Calvanese_EDFig1f_AGMHema.SeuratObj)
alias_table.G2M <- resolve_gene_aliases(g2m.genes, Calvanese_EDFig1f_AGMHema.SeuratObj)

match_map.S <- setNames(alias_table.S$matched, alias_table.S$original)
match_map.G2M <- setNames(alias_table.G2M$matched, alias_table.G2M$original)

s.genes_resolved <- unname(match_map.S)
g2m_resolved <- unname(match_map.G2M)

# Add module scores and phase assignments
Calvanese_EDFig1f_AGMHema.SeuratObj <- 
  CellCycleScoring(Calvanese_EDFig1f_AGMHema.SeuratObj, 
                   s.features = s.genes_resolved, 
                   g2m.features = g2m_resolved, 
                   set.ident = TRUE)

# View cell cycle scores and phase assignments
p <- DimPlot(Calvanese_EDFig1f_AGMHema.SeuratObj, group.by = "Phase",
             label = F, label.size = 7,
             cols = CC_COL) +
  ggtitle("UMAP with Cell Cycle Phase \nin Calvanese et al. AGM CS14&15 haematopoietic cells") +
  theme(axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1),
        axis.text = element_text(face = "bold", colour = "black"),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = .5, colour = "black", size = 15))

ggsave(
  filename = "plots/HM_AGM_CC_UMAP.pdf",
  p,
  width = 8, height = 7,
  device = cairo_pdf
)

# Histogram plot for CT_HM
df_cc <- Calvanese_EDFig1f_AGMHema.SeuratObj[[]] %>% 
  as.data.frame() %>% 
  dplyr::select(CT_HM, Phase) %>% 
  dplyr::filter(!is.na(CT_HM), !is.na(Phase)) %>%
  dplyr::group_by(CT_HM, Phase) %>% 
  dplyr::summarise(n = n(), .groups = "drop_last") %>% 
  dplyr::mutate(freq = n / sum(n)) %>% 
  dplyr::ungroup()

df_cc$CT_HM <- factor(df_cc$CT_HM, levels = CT_ORDER)

df_cc$Phase <- factor(df_cc$Phase, levels = c("G1", "S", "G2M"))

p <- ggplot(df_cc, aes(x = CT_HM, y = freq, fill = Phase)) +
  geom_col(color = "white", width = 0.8, linewidth = .1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Cell type", y = "Cell cycle phase (%)", fill = "Cell Cycle Phase") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major.x = element_blank()
  ) +
  scale_fill_manual(values = CC_COL) +
  ggtitle("Cell Cycle Phase Percentage \nin Calvanese et al. AGM CS14&15 haematopoietic cells") +
  theme(axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1), 
        axis.text = element_text(face = "bold", colour = "black"),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", hjust = .5, colour = "black", size = 17))


ggsave(
  filename = "plots/HM_AGM_CC_percentage.pdf",
  p,
  width = 9, height = 6,
  device = cairo_pdf
)

# Plot AGM & FL HSC together
library(readr)

## AGM HSC
df_agm_hsc <- Calvanese_EDFig1f_AGMHema.SeuratObj[[]] %>% 
  as.data.frame() %>% 
  dplyr::select(CT_HM, Phase) %>%
  dplyr::filter(CT_HM == "HSC", !is.na(Phase)) %>% 
  dplyr::group_by(Phase) %>% 
  dplyr::summarise(n = n(), .groups = "drop_last") %>% 
  dplyr::mutate(freq = n / sum(n),
                source = "AGM CS14-15 HSC") %>% 
  dplyr::ungroup()

## FL HSC
fl_cc <- read_csv("data/HSC_CCPhase_VCTime.csv")
## Get PCW5 barcode
HSC_bc <- read.csv("~/FL_wnn_cellmeta.v02.csv")
HSC_bc.v <- HSC_bc %>% 
  dplyr::filter(PCW == 5) %>% 
  dplyr::pull(X)

df_fl_hsc <- fl_cc %>% 
  dplyr::filter(barcode %in% HSC_bc.v) %>%
  dplyr::filter(!is.na(CC_Phase)) %>% 
  dplyr::mutate(
    CC_Phase = ifelse(CC_Phase == "G0", "G1", CC_Phase)
  ) %>% 
  dplyr::group_by(CC_Phase) %>% 
  dplyr::summarise(n = n(), .groups = "drop_last") %>% 
  dplyr::mutate(freq = n / sum(n),
                source = "FL PCW5 HSC") %>% 
  dplyr::ungroup() %>% 
  dplyr::rename(Phase = CC_Phase)

## Combine
df_plot <- bind_rows(df_agm_hsc, df_fl_hsc)

df_plot$Phase  <- factor(df_plot$Phase, levels = c("G1", "S", "G2M"))
df_plot$Phase <- droplevels(df_plot$Phase)
df_plot$source <- factor(df_plot$source, levels = c("AGM CS14-15 HSC", "FL PCW5 HSC"))

CC_COL <- c(
  "G1"   = "#2a6da1",
  "S"    = "#E6762B",
  "G2M"  = "#329848"
)

## Plot
p <- ggplot(df_plot, aes(x = source, y = freq, fill = Phase)) +
  geom_col(color = "white", width = 0.8, linewidth = .1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "", y = "Cell cycle phase (%)", 
       fill = "Cell Cycle Phase",
       title = "Cell cycle phase of AGM HSC vs FL HSC", 
       subtitle = "(AGM data is from Calvanese et al.)") +
  scale_fill_manual(values = CC_COL) +
  theme_classic() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title       = element_text(face = "bold", hjust = 0, size = 13),
        plot.subtitle    = element_text(face = "bold", hjust = 0, size = 9),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1))

ggsave(
  filename = "plots/AGM_vs_FL_HSC_CC_percentage.pdf",
  p,
  width = 5, height = 5,
  device = cairo_pdf
)


