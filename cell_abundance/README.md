# Cell Abundance and Developmental Changes in Cell Composition

This folder contains analyses of developmental changes in cell composition during human fetal liver hematopoiesis across different post-conceptional weeks (PCW).

## Overview

We employ two complementary approaches to identify and quantify cell type abundance changes across developmental stages:

1. **Milo (Spatial abstraction)** - Neighborhood-based differential abundance analysis
2. **GLMM (Generalized Linear Mixed Model)** - Sample-level differential abundance testing with statistical modeling

## Folder Structure

### `/milo`
Neighborhood-based differential abundance analysis using Milo.

**Key Scripts:**
- `01.prep_rna.R` - RNA data preprocessing, HVG selection, and Harmony batch correction
- `01.prep_atac.R` - ATAC data preprocessing with spectral decomposition and Harmony batch correction
- `02.prep_wnn.py` - Weighted Nearest Neighbor integration of RNA and ATAC modalities
- `03.milo_choose_k_p.py` - Parameter optimization for neighborhood size (k) and sampling proportion (p)
- `04.prep_milo_input.py` - Preparation of input data for Milo differential abundance testing
- `05.milo_test_da.py` - Execution of Milo differential abundance testing
- `06.annotate_nhoods.py` - Annotation of neighborhoods with cell type labels
- `07.milo_visualization.py` - Visualization of differential abundance results

**Outputs:**
- Neighborhood abundance estimates across PCW stages
- Cell type-annotated neighborhoods
- Statistical significance testing results
- Visualization plots (UMAP with neighborhood coloring, heatmaps, etc.)

### `/glmm`
Sample-level differential abundance analysis using Poisson generalized linear mixed models.

**Key Scripts:**
- `Ctype_DA_Poisson_bySample.Rmd` - Poisson GLMM model fitting with library/donor/batch as random effects and developmental stage/sex as fixed effects
- `Ctype_DA_Poisson_check.Rmd` - Model diagnostics and validation checks (residuals, overdispersion, etc.)
- `Ctype_DA_basic.Rmd` - Basic visualization of cell type composition across developmental stages
- `Ctype_DA_plot.Rmd` - Detailed plotting of differential abundance results with significance highlighting
- `CellTypeCompositionAnalysis.R` - Utility functions for count matrix creation and statistical modeling

**Outputs:**
- Cell type composition tables and plots
- Statistical test results (p-values, log fold changes)
- Model diagnostics and fit statistics
- Publication-quality figures showing developmental dynamics

## Methodological Approach

### Data Preparation
- RNA and ATAC sequencing data from fetal liver samples collected across PCW 8-24
- Quality control and filtering to remove low-quality cells
- Feature selection (HVGs for RNA, accessible peaks for ATAC)
- Batch correction using Harmony with library and donor as batch variables

### Statistical Testing

**Milo Analysis:**
- k-nearest neighbor graph construction (k=100) based on integrated RNA-ATAC space
- Non-overlapping neighborhood sampling
- Logistic regression differential abundance testing per neighborhood
- Cell type annotation of neighborhoods based on majority cell type

**GLMM Analysis:**
- Poisson generalized linear mixed model: `count ~ PCW * Sex + (1|Library) + (1|Donor) + (1|Batch)`
- PCW (developmental stage) and Sex as fixed effects of interest
- Library, Donor, and Batch as random effects to account for batch effects
- SHAP-based interpretation and RMT null distribution testing (when applicable)

### Key Covariates
- **PCW (Post-conceptional week):** Developmental stage (8-24 weeks)
- **Sex:** Male/Female biological sex
- **Batch:** Technical batch effects
- **Library:** RNA/ATAC library batch
- **Donor:** Individual donor effects

## Data Requirements

This analysis requires:
- Preprocessed and annotated single-cell/nucleus RNA-seq data (Seurat object format)
- Preprocessed and annotated single-cell ATAC-seq data (SnapATAC2 object format)
- Sample metadata including PCW, sex, donor ID, library ID, and batch information
- Cell type annotations (e.g., `anno_wnn_v51`)

## Output Interpretation

### Milo Results
- `logFC`: Log-fold change in neighborhood abundance between conditions
- `pval`: Raw p-value from logistic regression
- `FDR`: Corrected p-value (Benjamini-Hochberg)
- Neighborhoods with FDR < 0.05 indicate significant differential abundance

### GLMM Results
- **Main effects:** Developmental stage effects (PCW) and sex effects
- **Interactions:** PCW × Sex interactions indicate sex-specific developmental dynamics
- Statistical significance reported with p-values and confidence intervals
- Model assumptions checked via residual diagnostics

## References

For detailed methods, see the associated manuscript:
- **Temporal and regulatory landscape of human embryonic and foetal liver haematopoiesis mapped by single-nucleus multiomics**

### Software
- Milo: https://github.com/MarioniLab/miloR
- R packages: Seurat, lme4, tidyverse
- Python packages: anndata, scanpy, pertpy, muon, harmony

## Notes

- GLMM/SHAP analysis steps include stochastic components (XGBoost subsampling). Set the `seed` argument to reproduce exact results reported in the manuscript.
- All downstream steps (visualization, Jaccard analysis) are deterministic given fixed inputs.
- This analysis does not include automated clustering; cell types are defined based on prior annotation from the annotation pipeline.
