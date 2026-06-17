# HSC Subclustering

This folder contains the HSC-focused subclustering pipeline for fetal liver multiomic data. It includes RNA preprocessing, ATAC preprocessing, weighted nearest neighbor integration, multiVI modeling, cluster quality checks, and differential analysis support.

## Main scripts

- `01.prep_rna.R` — preprocess RNA counts, normalization, and prepare HSC RNA data

- `02.prep_atac.py` — preprocess ATAC counts, QC, dimension reduction

- `03.run_wnn.py` — integrate RNA and ATAC data with weighted nearest neighbors

- `04.run_multivi.py` — train multiVI on paired RNA/ATAC

- `05.basic_check.ipynb` — initial QC and data checks

- `06.choose_resolution.ipynb` — choose clustering resolution

- `07.ROGUE.R`: calculate ROGUE cluster purity metrics using RNA counts for several Leiden resolutions

- `08.plot_top_deg.ipynb` — plot top differential genes

- `09.composition_basic.Rmd` — composition analysis report

- `10.composition_milo/` prepares HSC data for Milo differential abundance testing after excluding samples with uncertain pcw information
  - `01.prep_atac.py` — ATAC preprocessing, Harmony batch correction, UMAP, and KNN.
  - `02.prep_wnn.py` — build WNN MuData, add HSC-specific subcluster labels, save `.h5mu`.
  - `03.milo_choose_k_p.py` — optimize Milo neighborhood parameters `k` and `p` via neighborhood size distributions.
  - `04.prep_milo_input.py` — export `adata_4_milo.h5ad` for Milo input.
  - `05.milo_wnn.ipynb` — run Milo test
  - `06.export_nhood_test_res.ipynb` — export the results as table

- `11.compsotion_Poisson.Rmd` — Poisson composition modeling

- `12.dev_de_lmm_subcluster.R` and `12.dev_de_lmm_subcluster.Rmd`: run linear mixed model differential analysis per subcluster.

