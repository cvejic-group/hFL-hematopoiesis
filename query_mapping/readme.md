# Query Mapping

This folder contains codes to do label transfer/query mapping from our PCW5 FL HSC to [public AGM data](https://www.nature.com/articles/s41586-022-04571-x) (`01.LabelTransfer_FL2Calvanese_PCW5HSC_AGM.R`) and characterize cell cycle phases in the same AGM data (`02.AGM_HSC_CC_HM.R`).

- `01.LabelTransfer_FL2Calvanese_PCW5HSC_AGM.R` - use the original cell type label and `symphone` package to map our PCW5 FL HSC onto CS14-15 AGM data.
- `02.AGM_HSC_CC_HM.R` - Using `Seurat` package to investigate the cell cycle phases of the same CS14-15 AGM data.