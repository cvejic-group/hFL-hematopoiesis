[![DOI](https://zenodo.org/badge/1081882728.svg)](https://doi.org/10.5281/zenodo.20024863)

# Multiomics analysis of fetal liver hematopoiesis

Here are scripts accompanying the manuscript: Temporal and regulatory landscape of human embryonic and foetal liver haematopoiesis mapped by single-nucleus multiomics

## Codebase

Code and scripts used for the custom analyses in the manuscript. These are shared to improve reproducibility of the main results, but the repository is not a complete, standalone workflow. Please see each sub-folder for more details.

* preprocess: data preprocessing and QC
  * 01.cellranger-arc
  * 02.souporcell
  * 03.scDblFinder
  * 04.qc
  * 05.infercnv

* annotation: cell type annotation
  * 01.wnn_anno
  * 02.multivi
  * 03.archr
  * 04.scavenge
    * preprocess
    * [scavenge-smk-pipeline](https://github.com/cvejic-group/scavenge-smk-devmult)
    * analysis

* cell_abundance: developmental changes in cell composition
  * milo
  * glmm

* trajectory: trajectory inference
  * moscot_temporal
  * multivi_diffusion
  * paga

* cell_cycle: cell-cycle related analysis
  * 01.qHSC
  * 02.Dynamic_CC_modeling

* peak_gene_links: identification of peak-gene links

* dev_differential: developmental effects on chromatin accessibility and gene expression
  * dev_diff_peak
  * dev_diff_gene

* hsc_subclustering: explore the heterogeneity within HSCs

* gwas_enrichment: GWAS signal enrichment using cell type-resolved regulatory elements
  * sldsc
  * magama

* DNAm: DNA methylation data analysis

* tf_footprinting: TF footprint identification with ChromBPNet, TF-MoDISCo, and Fi-NeMo
  * [chrombpnet-smk-pipeline](https://github.com/cvejic-group/chrombpnet-smk-devmult)
  * postprocess
  * analysis
  * [tobias-smk-pipeline](https://github.com/cvejic-group/tobias-smk-devmult)

* query_mapping: Mapping FL HSC to AGM data from [Calvanese et al. 2022, Nature](https://www.nature.com/articles/s41586-022-04571-x)

* grn_analysis: gene regulatory network (GRN) analysis with SCENIC+ and XGBoost model
  * scenicplus
  * celltype_resolved_grn
  * dev_diff_eRegulon: developmental effects on SCENIC+ eRegulons.
  * B_lineage_demonstration: A detailed demonstration of our eGRN analysis within B lineage.

* misc
  * variancePartition: variance partition analysis
  * cNMF: consensus non-negative matrix factorization
  * bam_frag_per_cell: generate BAM/fragment files for each cell type

* utils: small utilities, color palette



