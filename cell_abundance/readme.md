# About

This folder contains analyses of developmental changes in cell composition during human fetal liver hematopoiesis across different post-conceptional weeks (PCW).

## Overview

We used two complementary approaches to quantify cell type abundance changes across developmental stages:

1. **Milo (Spatial abstraction)** - Neighborhood-based differential abundance analysis
2. **GLMM (Generalized Linear Mixed Model)** - Sample-level differential abundance testing with statistical modeling

## Files

- `Ctype_DA_basic.Rmd` - Basic visualization of cell type composition across developmental stages
- `Ctype_DA_plot.Rmd` - Detailed plotting of differential abundance results with significance highlighting

### `/milo`

- `01.prep_rna.R` - RNA data preprocessing, HVG selection, and Harmony batch correction
- `01.prep_atac.R` - ATAC data preprocessing with spectral decomposition and Harmony batch correction
- `02.prep_wnn.py` - Weighted Nearest Neighbor integration of RNA and ATAC modalities
- `03.milo_choose_k_p.py` - Parameter optimization for neighborhood size (k) and sampling proportion (p)
- `04.prep_milo_input.py` - Preparation of input data for Milo differential abundance testing
- `05.milo_test_da.py` - Execution of Milo differential abundance testing
- `06.annotate_nhoods.py` - Annotation of neighborhoods with cell type labels
- `07.milo_visualization.py` - Visualization of differential abundance results

### `/glmm`

- `Ctype_DA_Poisson_check.Rmd` - Model diagnostics and validation checks (residuals, overdispersion, etc.)
- `Ctype_DA_Poisson_bySample.Rmd` - Poisson GLMM model fitting with library/donor/batch as random effects and developmental stage/sex as fixed effects
- `CellTypeCompositionAnalysis.R` - Utility functions for count matrix creation and statistical modeling

