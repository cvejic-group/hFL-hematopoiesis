## About

This folder collects multi-modal annotation, integration, and regulatory annotation workflows for fetal liver multiome data.

## Contents

- `01.wnn_anno`: RNA and ATAC preprocessing, joint WNN integration, iterative manual annotation, and cleaned joint dataset export.
- `02.multivi`: MultiVI model preparation, training, and visualization for joint RNA+ATAC integration across all PCWs.
- `03.archr`: ArchR-based ATAC annotation, peak calling, motif enrichment, peak-to-gene linking, and export utilities.
- `04.scavenge`: SCAVENGE GWAS trait relevance score workflows and ArchR-to-SummarizedExperiment export.


## Key compute environment

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

