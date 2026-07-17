#------- Initialization -------#

# Load essential packages
library(scater)
library(Seurat)
library(Signac)
library(ggplot2)
library(cowplot)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

# alias_resolve
library(org.Hs.eg.db)
library(AnnotationDbi)
resolve_gene_aliases <- function(genes, seurat_obj) {
  
  obj_genes <- rownames(seurat_obj)
  
  result <- lapply(genes, function(g) {
    
    # Pass 1: direct match
    if (g %in% obj_genes) {
      return(data.frame(original = g, matched = g, label = g, found = TRUE,
                        stringsAsFactors = FALSE))
    }
    
    # Pass 2: look up ENTREZID - try SYMBOL first, then ALIAS
    entrez <- suppressMessages(tryCatch(
      AnnotationDbi::mapIds(org.Hs.eg.db, keys = g,
                            column = "ENTREZID", keytype = "SYMBOL",
                            multiVals = "first"),
      error = function(e) NA_character_
    ))
    
    if (is.na(entrez)) {
      entrez <- suppressMessages(tryCatch(
        AnnotationDbi::mapIds(org.Hs.eg.db, keys = g,
                              column = "ENTREZID", keytype = "ALIAS",
                              multiVals = "first"),
        error = function(e) NA_character_
      ))
    }
    
    if (is.na(entrez)) {
      return(data.frame(original = g, matched = NA_character_, label = g,
                        found = FALSE, stringsAsFactors = FALSE))
    }
    
    # Get all aliases + official symbol for this gene
    all_aliases <- suppressMessages(tryCatch(
      AnnotationDbi::mapIds(org.Hs.eg.db, keys = entrez,
                            column = "ALIAS", keytype = "ENTREZID",
                            multiVals = "list")[[1]],
      error = function(e) character(0)
    ))
    official_symbol <- suppressMessages(tryCatch(
      AnnotationDbi::mapIds(org.Hs.eg.db, keys = entrez,
                            column = "SYMBOL", keytype = "ENTREZID",
                            multiVals = "first"),
      error = function(e) NA_character_
    ))
    
    candidates <- unique(c(official_symbol, all_aliases))
    candidates <- candidates[!is.na(candidates)]
    
    # Find which candidates exist in the Seurat object
    hits <- candidates[candidates %in% obj_genes]
    
    if (length(hits) == 0) {
      return(data.frame(original = g, matched = NA_character_, label = g,
                        found = FALSE, stringsAsFactors = FALSE))
    }
    
    matched <- hits[1]  # use first hit
    label   <- if (matched == g) g else paste0(g, "/", matched)
    
    return(data.frame(original = g, matched = matched, label = label,
                      found = TRUE, stringsAsFactors = FALSE))
  })
  
  do.call(rbind, result)
}


