# Cell Cycle

This folder contains HSC cell-cycle analysis split into two subfolders: quiescent HSC analysis (`01.qHSC`) and dynamic cell cycle modeling (`02.Dynamic_CC_modeling`).

## `01.qHSC`

- `01.CvsQ.DE_LMM.R` and `02.CvsQ.DE_LMM.Rmd` - use LMM to find genes differentially expressed between G0 and other HSCs.
- `03.CvsQ.DE_cmp_devDE.Rmd` - compares qHSC differential expression with developmental DE patterns.
- `04.CvsQ.DE_plot.Rmd` - visualizes qHSC differential expression with volcano plots and summary figures.
- `05.archr.HSC_CvsQ.chromVAR_diff_peak.R` and `05.archr.HSC_CvsQ.chromVAR_diff_peak.Rmd` - UseschromVAR to compare ATAC peak activity between G0 and other cells.
- `06.archr.HSC_CvsQ.chromVAR_ChIP_peak.R` and `06.archr.HSC_CvsQ.chromVAR_ChIP_peak.Rmd` - Extends chromVAR analysis to ChIP-derived peak sets.
- `07.monaLisa_metaplot.R` and `08.monaLisa_metaplot_ALLPeak.R` - Metaplots for promoter and distal regions between q-HSC and c-HSC. They are either from differential peaks, or from all peaks called in HSCs.
- `09.monaLisa_scatter.R` -  Scatter plot for TF motif enrichment analysis with monaLisa on q-HSC vs c-HSC grouped by promoter and distal.

## `02.Dynamic_CC_modeling`

- `CC_pseudotime_inference.ipynb` - Infers cell-cycle pseudotime for ordering cells along a cyclic trajectory.
- `CycSpline_fitting_functions.r` - Defines helper functions for cyclic spline and GAM modeling.
- `CycSpline_GEX_fitting.r` - Fits cyclic spline models to gene expression metacell data.
- `CycSpline_PEAK_fitting.r` - Fits cyclic spline models to accessibility peak metacell data.
- `Cyclic_GEX_clustering.r` - Clusters cyclic gene expression trajectories and creates pattern plots.
- `Cyclic_PEAK_annotation.r` - Annotates cyclic accessibility peaks with genomic context or regulatory labels.
- `Cyclic_PEAK_clustering.r` - Clusters cyclic peak patterns and generates heatmaps.
- `Cyclic_PEAK_tf_enrichment.r` - Performs TF motif enrichment on cyclic peak clusters.
- `Cyclic_PEAK_tf_PromoterVsDistalEnrichment.r` - Performs TF motif enrichment on cyclic peak grouped by promoter and distal.
- `Cyclic_PEAK_tf_PromoterVsDistalEnrichment_scatterplot.r` -  Scatter plot for TF motif enrichment on cyclic peak grouped by promoter and distal.

