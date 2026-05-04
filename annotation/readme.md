## About

This folder collects multi-modal annotation, integration, and regulatory annotation workflows for fetal liver multiome data.

## Contents

- `01.wnn_anno`: RNA and ATAC preprocessing, joint WNN integration, iterative manual annotation, and cleaned joint dataset export.
- `02.multivi`: MultiVI model preparation, training, and visualization for joint RNA+ATAC integration across all PCWs.
- `03.archr`: ArchR-based ATAC annotation, peak calling, motif enrichment, peak-to-gene linking, and export utilities.
- `04.scavenge`: SCAVENGE GWAS trait relevance score workflows and ArchR-to-SummarizedExperiment export.


## Key compute environment

### Seurat

```r
> sessionInfo()
R version 4.3.2 (2023-10-31)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 22.04.2 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
 [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

time zone: Europe/Berlin
tzcode source: system (glibc)

attached base packages:
[1] stats     graphics  grDevices datasets  utils     methods   base     

other attached packages:
[1] Seurat_5.0.1       SeuratObject_5.0.1 sp_2.1-2           workflowr_1.7.1   

loaded via a namespace (and not attached):
  [1] deldir_2.0-2           pbapply_1.7-2          gridExtra_2.3          rlang_1.1.6            magrittr_2.0.3        
  [6] git2r_0.33.0           RcppAnnoy_0.0.22       spatstat.geom_3.2-7    matrixStats_1.5.0      ggridges_0.5.4        
 [11] compiler_4.3.2         getPass_0.2-2          reshape2_1.4.4         png_0.1-8              callr_3.7.3           
 [16] vctrs_0.6.5            stringr_1.5.1          pkgconfig_2.0.3        fastmap_1.1.1          promises_1.3.2        
 [21] rmarkdown_2.29         ps_1.7.5               purrr_1.0.4            xfun_0.52              jsonlite_1.8.8        
 [26] goftest_1.2-3          later_1.3.2            spatstat.utils_3.1-0   irlba_2.3.5.1          parallel_4.3.2        
 [31] cluster_2.1.4          R6_2.5.1               ica_1.0-3              spatstat.data_3.0-3    stringi_1.8.2         
 [36] RColorBrewer_1.1-3     reticulate_1.34.0      parallelly_1.46.1      scattermore_1.2        lmtest_0.9-40         
 [41] Rcpp_1.0.14            knitr_1.50             tensor_1.5             future.apply_1.11.0    zoo_1.8-12            
 [46] sctransform_0.4.1      httpuv_1.6.13          Matrix_1.6-3           splines_4.3.2          igraph_1.5.1          
 [51] tidyselect_1.2.1       abind_1.4-8            rstudioapi_0.15.0      dichromat_2.0-0.1      yaml_2.3.7            
 [56] spatstat.random_3.2-2  spatstat.explore_3.2-5 codetools_0.2-20       miniUI_0.1.1.1         processx_3.8.2        
 [61] listenv_0.9.0          plyr_1.8.9             lattice_0.22-6         tibble_3.2.1           shiny_1.10.0          
 [66] ROCR_1.0-11            evaluate_1.0.4         Rtsne_0.17             future_1.69.0          fastDummies_1.7.3     
 [71] survival_3.5-7         polyclip_1.10-6        fitdistrplus_1.1-11    pillar_1.10.1          BiocManager_1.30.25   
 [76] whisker_0.4.1          KernSmooth_2.23-24     renv_1.0.3             plotly_4.10.3          generics_0.1.3        
 [81] rprojroot_2.0.4        RcppHNSW_0.6.0         ggplot2_3.5.1          scales_1.4.0           globals_0.19.0        
 [86] xtable_1.8-4           glue_1.6.2             lazyeval_0.2.2         tools_4.3.2            data.table_1.17.0     
 [91] RSpectra_0.16-2        RANN_2.6.1             fs_1.6.3               leiden_0.4.3.1         dotCall64_1.1-1       
 [96] cowplot_1.1.3          grid_4.3.2             tidyr_1.3.1            colorspace_2.1-1       nlme_3.1-163          
[101] patchwork_1.3.0.9000   cli_3.6.1              spatstat.sparse_3.0-3  spam_2.10-0            viridisLite_0.4.2     
[106] dplyr_1.1.4            uwot_0.2.3             gtable_0.3.6           digest_0.6.33          progressr_0.14.0      
[111] ggrepel_0.9.6          htmlwidgets_1.6.4      farver_2.1.2           htmltools_0.5.7        lifecycle_1.0.3       
[116] httr_1.4.7             here_1.0.1             mime_0.12              MASS_7.3-60.0.1
```

### MUON

```sh
$ session-info muon
muon	0.1.7
----	----
plotly	5.21.0
cffi	1.16.0
statsmodels	0.14.1
dill	0.3.8
tqdm	4.66.5
pytz	2024.1
llvmlite	0.43.0
python-dateutil	2.9.0.post0
texttable	1.7.0
numcodecs	0.12.1
patsy	0.5.6
msgpack	1.0.8
pycparser	2.22
pillow	10.3.0
zarr	2.17.2
igraph	0.11.4
wcwidth	0.2.5
umap-learn	0.5.6
joblib	1.4.0
cycler	0.12.1
natsort	8.4.0
PyYAML	6.0.1
scikit-learn	1.4.2
six	1.16.0
leidenalg	0.10.2
asciitree	0.3.3
pynndescent	0.5.12
kiwisolver	1.4.5
torch	2.2.2 (2.2.2+cu121)
psutil	5.9.0
pyarrow	16.0.0
numba	0.60.0
h5py	3.11.0
setuptools	68.2.2
matplotlib	3.8.4
----	----
Python	3.10.14 (main, Mar 21 2024, 16:24:04) [GCC 11.2.0]
OS	Linux-6.12.0-124.52.3.el10_1.x86_64-x86_64-with-glibc2.35
CPU	256/256 logical CPU cores, x86_64
Updated	2026-05-04 13:46
```

### scVI

```sh
$ session-info scvi
scvi-tools	1.2.2.post2
----	----
typing_extensions	4.12.2
zstandard	0.23.0
flax	0.10.2
jaraco.functools	4.0.1
lightning	2.5.0.post0
PyYAML	6.0.2
pillow	11.1.0
more-itertools	10.3.0
gmpy2	2.1.5
tqdm	4.67.1
colorama	0.4.6
torch	2.5.1.post207
pandas	2.2.3
jaraco.collections	5.1.0
numba	0.61.0
xarray	2025.1.2
scipy	1.15.2
charset-normalizer	3.4.1
kiwisolver	1.4.8
pytz	2024.1
lightning-utilities	0.12.0
psutil	5.9.0
jaraco.context	5.3.0
optax	0.2.4
msgpack	1.1.0
filelock	3.17.0
pyparsing	3.2.1
setuptools	75.8.0
Pygments	2.19.1
multipledispatch	0.6.0
natsort	8.4.0
numpy	2.1.3
fsspec	2025.2.0
session-info2	0.4.1
anndata	0.11.3
ml_collections	1.0.0
etils	1.12.0
mudata	0.3.1
packaging	24.2
toolz	1.0.0
ml-dtypes	0.5.1
h5py	3.13.0
cycler	0.12.1
jax-cuda12-plugin	0.4.35
numpyro	0.17.0
jax	0.4.35
jax-cuda12-pjrt	0.4.35
jaxlib	0.4.35
llvmlite	0.44.0
matplotlib	3.10.0
six	1.17.0
chex	0.1.88
setuptools-scm	8.1.0
pyro-ppl	1.9.1+ab0491a
docrep	0.3.2
sparse	0.15.5
rich	13.9.4
opt_einsum	3.4.0
absl-py	2.1.0
mpmath	1.3.0
scikit-learn	1.5.2
python-dateutil	2.9.0.post0
torchmetrics	1.6.1
joblib	1.4.2
threadpoolctl	3.5.0
jaraco.text	3.12.1
numexpr	2.10.1
sympy	1.13.3
----	----
Python	3.12.9 | packaged by conda-forge | (main, Feb 14 2025, 08:00:06) [GCC 13.3.0]
OS	Linux-6.12.0-124.52.3.el10_1.x86_64-x86_64-with-glibc2.35
CPU	256/256 logical CPU cores, x86_64
Updated	2026-05-04 13:52
```

### ArchR

```r
> sessionInfo()
R version 4.1.2 (2021-11-01)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04.3 LTS

Matrix products: default
BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.8.so

Random number generation:
 RNG:     L'Ecuyer-CMRG 
 Normal:  Inversion 
 Sample:  Rejection 
 
locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
 [1] parallel  stats4    grid      stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] BSgenome.Hsapiens.UCSC.hg38_1.4.4 BSgenome_1.62.0                   rtracklayer_1.54.0               
 [4] Biostrings_2.62.0                 XVector_0.34.0                    rhdf5_2.38.1                     
 [7] SummarizedExperiment_1.24.0       Biobase_2.54.0                    RcppArmadillo_0.12.8.4.0         
[10] Rcpp_1.0.12                       Matrix_1.5-0                      GenomicRanges_1.46.1             
[13] GenomeInfoDb_1.30.1               IRanges_2.28.0                    S4Vectors_0.32.4                 
[16] BiocGenerics_0.40.0               sparseMatrixStats_1.6.0           MatrixGenerics_1.6.0             
[19] matrixStats_1.1.0                 data.table_1.15.4                 stringr_1.5.1                    
[22] plyr_1.8.9                        magrittr_2.0.3                    ggplot2_3.4.0                    
[25] gtable_0.3.5                      gtools_3.9.5                      gridExtra_2.3                    
[28] devtools_2.4.5                    usethis_2.2.3                     ArchR_1.0.3                      
[31] readr_2.1.5                       dplyr_1.1.4                       workflowr_1.7.1                  

loaded via a namespace (and not attached):
 [1] colorspace_2.1-0         rjson_0.2.21             ellipsis_0.3.2           rsconnect_0.8.25         rprojroot_2.0.4         
 [6] fs_1.6.4                 rstudioapi_0.16.0        remotes_2.5.0            bit64_4.0.5              fansi_1.0.6             
[11] cachem_1.1.0             knitr_1.48               pkgload_1.4.0            Rsamtools_2.10.0         Cairo_1.6-2             
[16] shiny_1.8.1.1            compiler_4.1.2           httr_1.4.7               fastmap_1.2.0            cli_3.6.3               
[21] later_1.3.2              htmltools_0.5.8.1        tools_4.1.2              glue_1.7.0               GenomeInfoDbData_1.2.7  
[26] vctrs_0.6.5              rhdf5filters_1.6.0       xfun_0.45                ps_1.7.7                 mime_0.12               
[31] miniUI_0.1.1.1           lifecycle_1.0.4          restfulr_0.0.15          XML_3.99-0.17            getPass_0.2-4           
[36] zlibbioc_1.40.0          scales_1.3.0             vroom_1.6.5              hms_1.1.3                promises_1.3.0          
[41] yaml_2.3.9               memoise_2.0.1            stringi_1.8.4            BiocIO_1.4.0             pkgbuild_1.4.4          
[46] BiocParallel_1.28.3      rlang_1.1.4              pkgconfig_2.0.3          bitops_1.0-7             evaluate_0.24.0         
[51] lattice_0.22-6           purrr_1.0.2              Rhdf5lib_1.16.0          GenomicAlignments_1.30.0 htmlwidgets_1.6.4       
[56] bit_4.0.5                processx_3.8.4           tidyselect_1.2.1         here_1.0.1               R6_2.5.1                
[61] generics_0.1.3           profvis_0.3.8            DelayedArray_0.20.0      pillar_1.9.0             whisker_0.4.1           
[66] withr_3.0.0              RCurl_1.98-1.14          tibble_3.2.1             crayon_1.5.3             utf8_1.2.4              
[71] tzdb_0.4.0               rmarkdown_2.27           urlchecker_1.0.1         callr_3.7.6              git2r_0.33.0            
[76] digest_0.6.36            xtable_1.8-4             httpuv_1.6.15            munsell_0.5.1            sessioninfo_1.2.2
```

