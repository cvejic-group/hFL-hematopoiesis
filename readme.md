
# Multiomics analysis of fetal liver hematopoiesis

Here are scripts accompanying the manuscript: Temporal and regulatory landscape of human foetal liver haematopoiesis mapped by single-nucleus multiomics


## Contents

Code and scripts used for the custom analyses in the manuscript. These are shared to improve reproducibility of the main results, but the repository is not a complete, standalone workflow.

* preprocess: data processing
  * 01.cellranger-arc
  * 02.souporcell
  * 03.doublets
  * 04.qc
  * 05.infercnv

* annotation: cell type annotation
  * 01.wnn_anno
  * 02.multivi
  * 03.archr

* cell_abundance: developmental changes in cell composition
  * milo
  * glmm

* trajectory: trajectory inference
  * moscot_temporal
  * multivi_diffusion
  * paga

* cell_cycle: cell-cycle related analysis
  * 01.qHSC

* peak_gene_links: identification of peak-gene links

* dev_differential: developmental effects on chromatin accessibility and gene expression
  * dev_diff_peak
  * dev_diff_gene

* gwas_enrichment: GWAS signal enrichment using cell type-resolved regulatory elements
  * sldsc
  * magama

* DNAm: DNA methylation

* misc
  * variancePartition
  * cNMF
  * bam_frag_per_cell

* utils


