# About

This folder contains code for developmental differential analysis of HSCs using linear mixed models (LMM) on both chromatin accessibility and gene expression.

The LMM method was first described by [Young et al. 2021](https://www.nature.com/articles/s41588-021-00875-2) and also employed by other studies. Statistical significance was estimated by local true sign rate (LTSR). Please find more details in the Supplementary Notes of [Young et al. 2021](https://www.nature.com/articles/s41588-021-00875-2).

## Contents

- `LMM.R`: core LMM implementation adapted from the [SKM_ageing_atlas](https://github.com/Teichlab/SKM_ageing_atlas) repository.
  - `LICENSE`: a copy of the LICENSE file from the [SKM_ageing_atlas](https://github.com/Teichlab/SKM_ageing_atlas) repository.

- `HSC.prep_temporal.R`: prepares input matrices and metadata for temporal analysis, including metacell aggregation and Seurat normalization/scaling.

- `dev_diff_gene`: gene-level developmental differential expression analysis.
  - `HSC.dev_de_lmm.R`: core script
  - `HSC.dev_de_lmm.Rmd`: functional enrichment and visualization
  - `HSC.dev_de_plot.Rmd`: visualization
  - `HSC.dev_de_lmm_qScoreAdj.r`: adjust for quiescence
  - `HSC.dev_de_lmm_qScoreAdj.rmd`: comparison with model without adjusting for quiescence

- `dev_diff_peak`: peak-level differential accessibility analysis, metacell-level modeling, and TF motif enrichment.
  - `HSC.dev_da_metacell.R`: core script
  - `HSC.dev_da_metacell.Rmd`: peak annotation, functional enrichment and visualization
  - `HSC.dev_da_CpG.Rmd`: GC content check of differential peaks
  - `HSC.dev_da_monaLisa.Rmd`: TF motif enrichment
  - `HSC.dev_da_plot.Rmd`: visualization

## Notes

- The analysis workflow is based on LMMs for both chromatin accessibility and gene expression, with explicit control for effects of library, donor, sex, and batch etc.

- Peak-level modeling is performed on **metacells** to reduce sparsity and improve statistical power.

- Motif enrichment is derived from binned peak effect sizes and filtered for strong significance and biological support.

