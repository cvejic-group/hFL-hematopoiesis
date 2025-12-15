################################################################
### Check eGRN Filtering Results with GO Enrichment Analysis ###
################################################################

###-------------------Note-------------------###
###-------      Run with R-4.4.2      -------###
###------------------------------------------###

##################
# Initialization #
##################

setwd("~/work/")
source("./00.Initialization.R")

# Load additional packages
library(AUCell)
library(ArchR)
library(clusterProfiler)
library(org.Hs.eg.db)

###################################################
### Load Cell Type specific eGRN by XGBoost RMT ###
###################################################

# Load CT-eRegulons info
CT_eReg_XGBoost_RMT.df <- readRDS(paste0(RES_DIR, "XGBoost_RMT_TFRes.rds"))

# Load SCENIC+ eRegulon metadata
raw_eRegulon_meta.df <- readRDS(paste0(DATA_DIR, "eRegulon_raw_meta.rds"))
BA_eRegulon_meta.df <- raw_eRegulon_meta.df %>% dplyr::filter(class == "+/+")

##################################
### Do GO for CT-eGRN TF + TGs ###
##################################

# Get CT gene list
CT_eGRN.l <- list()
for (CT in CT_ORDER) {
  CT_eGRN.l[[CT]] <- BA_eRegulon_meta.df %>% 
    dplyr::filter(TF %in% CT_eRegulon.l[[CT]]) %>% 
    pull(Gene) %>% 
    unique %>% 
    union(CT_eRegulon.l[[CT]], .)
}

# Do GO
## Covert SYMBOL to ENTREZID
CT_entrez.l <- lapply(CT_eGRN.l, function(sym) {
  df <- bitr(sym,
             fromType = "SYMBOL",
             toType = "ENTREZID",
             OrgDb = org.Hs.eg.db)
  unique(df$ENTREZID)
})

## GO Biological Process enrichment
go_raw.l <- lapply(CT_entrez.l, function(eids) {
  enrichGO(gene = eids,
           OrgDb = org.Hs.eg.db,
           keyType = "ENTREZID",
           ont = "BP",
           pAdjustMethod = "BH",
           qvalueCutoff  = 0.05)
})

# ## Use simplify to remove redundant terms
# go_simp.l <- lapply(go_raw.l, function(x) {
#   if (is.null(x) || nrow(x@result)==0) return(x)
#   simplify(x,
#            cutoff  = 0.7,
#            by = "p.adjust",
#            select_fun = min)
# })

## Generate dotplots
for (ct in names(go_raw.l)) {
  go.obj <- go_raw.l[[ct]]
  if (is.null(go.obj) || nrow(go.obj@result)==0) {
    message("No significant GO for ", ct)
    next
  }
  p <- dotplot(go.obj, showCategory = 20) +
    ggtitle(paste0(ct, " GO BP")) +
    theme_minimal() +
    theme(
      plot.title = element_text(face="bold", hjust = 0.5)
    )
  ggsave(filename = paste0(FIG_DIR, "GO_BP_", ct, ".pdf"), 
         plot = p, width = 6, height = 10)
}



