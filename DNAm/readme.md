# About

This folder contains the DNA methylation analysis pipeline for developmental methylation profiling and integration with chromatin and gene expression features.

The pipeline starts from Nanopore methylation calls generated with Dorado and the [Epi2me human variation workflow (v2.7.2)](https://github.com/epi2me-labs/wf-human-variation/tree/v2.7.2), producing `bedmethyl` files from `modBAM` output.

## Contents

- `01.load_mat_from_bed.R`: imports 5mC data from Nanopore-generated `bedmethyl.gz` files into `bsseq` objects.
- `02.filter_and_smooth.R`: filters CpG sites with missing coverage, smooths methylation values using BSmooth, and retains autosomal CpGs for downstream analysis.
- `03.basic_check_density.Rmd`: checks methylation density distributions and quality control.
- `03.basic_check_pca.Rmd`: performs PCA on smoothed methylation data and visualizes donor-level sample structure.
- `04.genomic_elements.R`: computes regional methylation summaries for promoters, ATAC peaks, ENCODE cCREs, and EpiMap K562 chromatin state segments.
- `05.get_atac_pb_by_donor.R`: generates donor-level pseudobulk ATAC peak matrices from ArchR.
- `05.get_rna_pb_by_donor.R`: generates donor-level pseudobulk gene expression matrices from Seurat.
- `06.DSS.R` / `06.DSS.Rmd`: fits DSS models and calls DMRs while accounting for PCW, sex, and Nanopore batch.
- `07.ChromVAR_DMR.R` / `07.ChromVAR_DMR.Rmd`: evaluates ChromVAR accessibility scores for hyper- and hypomethylated DMRs in FL HSCs.
- `08.DMR.olp_HSC_devDA_peaks.Rmd`: integrates DMRs with HSC differential accessibility peaks and performs overlap analyses.
- `DNAm.plot.Rmd`: creates summary visualizations of methylation PCA, chromatin state relationships, and DMR enrichments.

