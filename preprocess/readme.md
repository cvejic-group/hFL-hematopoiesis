
This folder contains the preprocessing pipeline used for snMultiome data demultiplexing, quality control, doublet detection, and *in silico* karyotype inference.

The pipeline is organized in the following subfolders:

- `01.cellranger-arc`: run Cell Ranger ARC alignment and generate gene expression and ATAC peak matrices.
- `02.souporcell`: genotype-based demultiplexing for pooled libraries using souporcell.
- `03.scDblFinder`: RNA doublet detection using `scDblFinder`.
- `04.qc`: QC metric calculation and joint cell filtering across RNA and ATAC modalities.
- `05.infercnv`: inferCNV-based karyotype inference for samples without tissue-bank karyotype data.

