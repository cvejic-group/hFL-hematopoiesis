# SCAVENGE Processing
This folder contains processing and analysis scripts for exploring GWAS trait-associated SNP enrichment in the atlas cell types.

## Input Files
The initial [`preprocess/01-export-archr-to-se-3.Rmd`](04.scavenge/preprocess/01-export-archr-to-se-3.Rmd) script begins from the ArchR project files, after generation of the post-Harmony LSI (ATAC) embedding,
which is computed in [`03.archr/archr.02.addUMAP.R`](03.archr/archr.02.addUMAP.R). 
Briefly, the LSI embedding was generated using `addIterativeLSI(...)` with 2 iterations, 25K variable features, 30 dimensions, and a `corCutOff` of 0.5. 
This was then followed by `addHarmony(...)` on this layer, using the same `corCutOff` and grouping by `libraryID` and `donorID`.

The resulting `SummarizedExperiment` object is then used in the **cvejic-group/scavenge-smk-devmult** pipeline to compute the SCAVENGE metrics (TRS, enrichment).

The `SummarizedExperiment` object and the per cell TRS results are used in the [`analysis/sfig3b-scavenge-baseline-trait-set.Rmd`](04.scavenge/analysis/sfig3b-scavenge-baseline-trait-set.Rmd) script to the summary figure.
The script specifically used the pipeline output `results/trs/df_trs_scavenge_yu_2022_knn30_bg200_rng20241106.csv.gz`.

## Software Environments
The following software environments were used:

### Preprocessing (ArchR)
**`sessionInfo()` output**

<details>

```r
## R version 4.1.3 (2022-03-10)
## Platform: x86_64-conda-linux-gnu (64-bit)
## Running under: Ubuntu 24.04 LTS
## 
## Matrix products: default
## BLAS/LAPACK: /work/aaa/miniforge3/envs/archr-dev-r41/lib/libopenblasp-r0.3.27.so
## 
## Random number generation:
##  RNG:     L'Ecuyer-CMRG 
##  Normal:  Inversion 
##  Sample:  Rejection 
##  
## locale:
##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
## 
## attached base packages:
##  [1] parallel  grid      stats4    stats     graphics  grDevices utils    
##  [8] datasets  methods   base     
## 
## other attached packages:
##  [1] BiocParallel_1.28.3               BSgenome.Hsapiens.UCSC.hg38_1.4.4
##  [3] BSgenome_1.62.0                   rtracklayer_1.54.0               
##  [5] Biostrings_2.62.0                 XVector_0.34.0                   
##  [7] RcppArmadillo_0.12.4.0.0          Rcpp_1.0.10                      
##  [9] sparseMatrixStats_1.6.0           data.table_1.14.8                
## [11] stringr_1.5.0                     plyr_1.8.8                       
## [13] ggplot2_3.4.2                     gtable_0.3.3                     
## [15] gtools_3.9.4                      gridExtra_2.3                    
## [17] devtools_2.4.5                    usethis_2.2.0                    
## [19] ArchR_1.0.3                       SummarizedExperiment_1.24.0      
## [21] Biobase_2.54.0                    GenomicRanges_1.46.1             
## [23] GenomeInfoDb_1.30.1               HDF5Array_1.22.1                 
## [25] rhdf5_2.38.1                      DelayedArray_0.20.0              
## [27] IRanges_2.28.0                    S4Vectors_0.32.4                 
## [29] MatrixGenerics_1.6.0              matrixStats_1.0.0                
## [31] BiocGenerics_0.40.0               Matrix_1.5-4.1                   
## [33] magrittr_2.0.3                   
## 
## loaded via a namespace (and not attached):
##  [1] bitops_1.0-7             fs_1.6.2                 tools_4.1.3             
##  [4] profvis_0.3.8            bslib_0.5.0              utf8_1.2.3              
##  [7] R6_2.5.1                 colorspace_2.1-0         rhdf5filters_1.6.0      
## [10] urlchecker_1.0.1         withr_2.5.0              tidyselect_1.2.0        
## [13] prettyunits_1.1.1        processx_3.8.1           compiler_4.1.3          
## [16] cli_3.6.1                sass_0.4.6               scales_1.2.1            
## [19] callr_3.7.3              Rsamtools_2.10.0         digest_0.6.31           
## [22] rmarkdown_2.22           pkgconfig_2.0.3          htmltools_0.5.5         
## [25] sessioninfo_1.2.2        fastmap_1.1.1            htmlwidgets_1.6.2       
## [28] rlang_1.1.1              shiny_1.7.4              BiocIO_1.4.0            
## [31] jquerylib_0.1.4          generics_0.1.3           jsonlite_1.8.5          
## [34] dplyr_1.1.2              RCurl_1.98-1.12          GenomeInfoDbData_1.2.7  
## [37] munsell_0.5.0            Rhdf5lib_1.16.0          fansi_1.0.4             
## [40] lifecycle_1.0.3          stringi_1.7.12           yaml_2.3.7              
## [43] zlibbioc_1.40.0          pkgbuild_1.4.0           promises_1.2.0.1        
## [46] crayon_1.5.2             miniUI_0.1.1.1           lattice_0.21-8          
## [49] knitr_1.43               ps_1.7.5                 pillar_1.9.0            
## [52] rjson_0.2.21             pkgload_1.3.2            XML_3.99-0.14           
## [55] glue_1.6.2               evaluate_0.21            remotes_2.4.2           
## [58] vctrs_0.6.2              httpuv_1.6.11            purrr_1.0.1             
## [61] cachem_1.0.8             xfun_0.39                mime_0.12               
## [64] xtable_1.8-4             restfulr_0.0.15          later_1.3.1             
## [67] tibble_3.2.1             GenomicAlignments_1.30.0 memoise_2.0.1           
## [70] ellipsis_0.3.2
```
  
</details>

**Conda Environment**
<details>

```yaml
name: archr-dev-r41
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - _libgcc_mutex=0.1=conda_forge
  - _openmp_mutex=4.5=2_gnu
  - _r-mutex=1.0.1=anacondar_1
  - binutils=2.43=h4852527_1
  - binutils_impl_linux-64=2.43=h4bf12b8_1
  - binutils_linux-64=2.43=h4852527_1
  - bioconductor-annotate=1.72.0=r41hdfd78af_0
  - bioconductor-annotationdbi=1.56.2=r41hdfd78af_0
  - bioconductor-biobase=2.54.0=r41hc0cfd56_2
  - bioconductor-biocgenerics=0.40.0=r41hdfd78af_0
  - bioconductor-biocio=1.4.0=r41hdfd78af_0
  - bioconductor-biocparallel=1.28.3=r41hc247a5b_1
  - bioconductor-biostrings=2.62.0=r41hc0cfd56_2
  - bioconductor-bsgenome=1.62.0=r41hdfd78af_0
  - bioconductor-bsgenome.hsapiens.ucsc.hg38=1.4.4=r41hdfd78af_1
  - bioconductor-chromvar=1.16.0=r41hc247a5b_2
  - bioconductor-cner=1.30.0=r41hc0cfd56_2
  - bioconductor-complexheatmap=2.10.0=r41hdfd78af_0
  - bioconductor-delayedarray=0.20.0=r41hc0cfd56_2
  - bioconductor-dirichletmultinomial=1.36.0=r41hda872b5_3
  - bioconductor-genomeinfodb=1.30.1=r41hdfd78af_0
  - bioconductor-genomeinfodbdata=1.2.7=r41hdfd78af_2
  - bioconductor-genomicalignments=1.30.0=r41hc0cfd56_2
  - bioconductor-genomicranges=1.46.1=r41hc0cfd56_1
  - bioconductor-go.db=3.14.0=r41hdfd78af_1
  - bioconductor-hdf5array=1.22.1=r41hc0cfd56_1
  - bioconductor-iranges=2.28.0=r41hc0cfd56_2
  - bioconductor-keggrest=1.34.0=r41hdfd78af_0
  - bioconductor-matrixgenerics=1.6.0=r41hdfd78af_0
  - bioconductor-motifmatchr=1.16.0=r41hc247a5b_2
  - bioconductor-rhdf5=2.38.1=r41hbe1951d_0
  - bioconductor-rhdf5filters=1.6.0=r41hc247a5b_2
  - bioconductor-rhdf5lib=1.16.0=r41hc0cfd56_2
  - bioconductor-rhtslib=1.26.0=r41hc0cfd56_2
  - bioconductor-rsamtools=2.10.0=r41hc247a5b_2
  - bioconductor-rtracklayer=1.54.0=r41h171f361_4
  - bioconductor-s4vectors=0.32.4=r41hc0cfd56_0
  - bioconductor-seqlogo=1.60.0=r41hdfd78af_0
  - bioconductor-sparsematrixstats=1.6.0=r41hc247a5b_2
  - bioconductor-summarizedexperiment=1.24.0=r41hdfd78af_0
  - bioconductor-tfbstools=1.32.0=r41hc0cfd56_2
  - bioconductor-xvector=0.34.0=r41hc0cfd56_2
  - bioconductor-zlibbioc=1.40.0=r41hc0cfd56_2
  - bwidget=1.9.14=ha770c72_1
  - bzip2=1.0.8=h4bc722e_7
  - c-ares=1.34.2=heb4867d_0
  - c-compiler=1.8.0=h2b85faf_1
  - ca-certificates=2024.8.30=hbcca054_0
  - cairo=1.16.0=ha61ee94_1014
  - compilers=1.8.0=ha770c72_1
  - curl=7.86.0=h7bff187_1
  - cxx-compiler=1.8.0=h1a2810e_1
  - expat=2.6.3=h5888daf_0
  - font-ttf-dejavu-sans-mono=2.37=hab24e00_0
  - font-ttf-inconsolata=3.000=h77eed37_0
  - font-ttf-source-code-pro=2.038=h77eed37_0
  - font-ttf-ubuntu=0.83=h77eed37_3
  - fontconfig=2.14.2=h14ed4e7_0
  - fonts-conda-ecosystem=1=0
  - fonts-conda-forge=1=0
  - fortran-compiler=1.8.0=h36df796_1
  - freetype=2.12.1=h267a509_2
  - fribidi=1.0.10=h36c2ea0_0
  - gcc=13.3.0=h9576a4e_1
  - gcc_impl_linux-64=13.3.0=hfea6d02_1
  - gcc_linux-64=13.3.0=hc28eda2_4
  - geos=3.11.2=hcb278e6_0
  - gettext=0.22.5=he02047a_3
  - gettext-tools=0.22.5=he02047a_3
  - gfortran=13.3.0=h9576a4e_1
  - gfortran_impl_linux-64=13.3.0=h10434e7_1
  - gfortran_linux-64=13.3.0=hb919d3a_4
  - glpk=5.0=h445213a_0
  - gmp=6.3.0=hac33072_2
  - graphite2=1.3.13=h59595ed_1003
  - gsl=2.7=he838d99_0
  - gxx=13.3.0=h9576a4e_1
  - gxx_impl_linux-64=13.3.0=hdbfa832_1
  - gxx_linux-64=13.3.0=h6834431_4
  - harfbuzz=6.0.0=h8e241bc_0
  - icu=70.1=h27087fc_0
  - jpeg=9e=h0b41bf4_3
  - kernel-headers_linux-64=3.10.0=he073ed8_17
  - keyutils=1.6.1=h166bdaf_0
  - krb5=1.19.3=h3790be6_0
  - ld_impl_linux-64=2.43=h712a8e2_1
  - lerc=4.0.0=h27087fc_0
  - libasprintf=0.22.5=he8f35ee_3
  - libasprintf-devel=0.22.5=he8f35ee_3
  - libblas=3.9.0=24_linux64_openblas
  - libcblas=3.9.0=24_linux64_openblas
  - libcurl=7.86.0=h7bff187_1
  - libdeflate=1.17=h0b41bf4_0
  - libedit=3.1.20191231=he28a2e2_2
  - libev=4.33=hd590300_2
  - libexpat=2.6.3=h5888daf_0
  - libffi=3.4.2=h7f98852_5
  - libgcc=14.2.0=h77fa898_1
  - libgcc-devel_linux-64=13.3.0=h84ea5a7_101
  - libgcc-ng=14.2.0=h69a702a_1
  - libgettextpo=0.22.5=he02047a_3
  - libgettextpo-devel=0.22.5=he02047a_3
  - libgfortran=14.2.0=h69a702a_1
  - libgfortran-ng=14.2.0=h69a702a_1
  - libgfortran5=14.2.0=hd5240d6_1
  - libgit2=1.5.1=ha98c156_0
  - libglib=2.78.1=hebfc3b9_0
  - libgomp=14.2.0=h77fa898_1
  - libiconv=1.17=hd590300_2
  - liblapack=3.9.0=24_linux64_openblas
  - libnghttp2=1.51.0=hdcd2b5c_0
  - libopenblas=0.3.27=pthreads_hac2b453_1
  - libpng=1.6.43=h2797004_0
  - libsanitizer=13.3.0=heb74ff8_1
  - libsodium=1.0.18=h36c2ea0_1
  - libssh2=1.10.0=haa6b8db_3
  - libstdcxx=14.2.0=hc0a3c3a_1
  - libstdcxx-devel_linux-64=13.3.0=h84ea5a7_101
  - libstdcxx-ng=14.2.0=h4852527_1
  - libtiff=4.5.0=h6adf6a1_2
  - libuuid=2.38.1=h0b41bf4_0
  - libwebp-base=1.4.0=hd590300_0
  - libxcb=1.13=h7f98852_1004
  - libxml2=2.10.3=hca2bb57_4
  - libzlib=1.2.13=h4ab18f5_6
  - make=4.4.1=hb9d3cd8_2
  - ncurses=6.5=he02047a_1
  - openssl=1.1.1w=hd590300_0
  - pandoc=2.19.2=h32600fe_2
  - pango=1.50.14=hd33c08f_0
  - pcre2=10.40=hc3806b6_0
  - pixman=0.43.2=h59595ed_0
  - pthread-stubs=0.4=hb9d3cd8_1002
  - r-abind=1.4_5=r41hc72bb7e_1004
  - r-askpass=1.1=r41h06615bd_3
  - r-assertthat=0.2.1=r41hc72bb7e_3
  - r-backports=1.4.1=r41h06615bd_1
  - r-base=4.1.3=h2f963a2_5
  - r-base64enc=0.1_3=r41h06615bd_1005
  - r-bh=1.81.0_1=r41hc72bb7e_0
  - r-bit=4.0.5=r41h06615bd_0
  - r-bit64=4.0.5=r41h06615bd_1
  - r-bitops=1.0_7=r41h06615bd_1
  - r-blob=1.2.4=r41hc72bb7e_0
  - r-brew=1.0_8=r41hc72bb7e_1
  - r-brio=1.1.3=r41h06615bd_1
  - r-broom=1.0.5=r41hc72bb7e_0
  - r-bslib=0.5.0=r41hc72bb7e_0
  - r-cachem=1.0.8=r41h57805ef_0
  - r-callr=3.7.3=r41hc72bb7e_0
  - r-catools=1.18.2=r41h7525677_1
  - r-cellranger=1.1.0=r41hc72bb7e_1005
  - r-chromvarmotifs=0.2.0=r41hdfd78af_0
  - r-circlize=0.4.15=r41hc72bb7e_1
  - r-cli=3.6.1=r41h38f115c_0
  - r-clipr=0.8.0=r41hc72bb7e_1
  - r-clue=0.3_64=r41h133d619_0
  - r-cluster=2.1.4=r41h8da6f51_0
  - r-codetools=0.2_19=r41hc72bb7e_0
  - r-colorspace=2.1_0=r41h133d619_0
  - r-commonmark=1.9.0=r41h133d619_0
  - r-conflicted=1.2.0=r41h785f33e_0
  - r-cowplot=1.1.1=r41hc72bb7e_1
  - r-cpp11=0.4.7=r41hc72bb7e_0
  - r-crayon=1.5.2=r41hc72bb7e_1
  - r-credentials=1.3.2=r41hc72bb7e_1
  - r-crosstalk=1.2.0=r41hc72bb7e_1
  - r-curl=4.3.3=r41h06615bd_1
  - r-data.table=1.14.8=r41h133d619_0
  - r-dbi=1.1.3=r41hc72bb7e_1
  - r-dbplyr=2.3.2=r41hc72bb7e_0
  - r-deldir=1.0_9=r41h61816a4_0
  - r-desc=1.4.2=r41hc72bb7e_1
  - r-devtools=2.4.5=r41hc72bb7e_1
  - r-diffobj=0.3.5=r41h06615bd_1
  - r-digest=0.6.31=r41h38f115c_0
  - r-doparallel=1.0.17=r41hc72bb7e_1
  - r-dotcall64=1.0_2=r41hac0b197_1
  - r-downlit=0.4.2=r41hc72bb7e_1
  - r-dplyr=1.1.2=r41ha503ecb_0
  - r-dqrng=0.3.0=r41h7525677_1
  - r-dt=0.28=r41hc72bb7e_0
  - r-dtplyr=1.3.1=r41hc72bb7e_0
  - r-ellipsis=0.3.2=r41h06615bd_1
  - r-evaluate=0.21=r41hc72bb7e_0
  - r-fansi=1.0.4=r41h133d619_0
  - r-farver=2.1.1=r41h7525677_1
  - r-fastdummies=1.6.3=r41hc72bb7e_1
  - r-fastmap=1.1.1=r41h38f115c_0
  - r-fitdistrplus=1.1_11=r41hc72bb7e_0
  - r-fnn=1.1.3.2=r41h38f115c_0
  - r-fontawesome=0.5.1=r41hc72bb7e_0
  - r-forcats=1.0.0=r41hc72bb7e_0
  - r-foreach=1.5.2=r41hc72bb7e_1
  - r-formatr=1.14=r41hc72bb7e_0
  - r-fs=1.6.2=r41ha503ecb_0
  - r-futile.logger=1.4.3=r41hc72bb7e_1004
  - r-futile.options=1.0.1=r41hc72bb7e_1003
  - r-future=1.32.0=r41hc72bb7e_0
  - r-future.apply=1.11.0=r41hc72bb7e_0
  - r-gargle=1.4.0=r41h785f33e_0
  - r-generics=0.1.3=r41hc72bb7e_1
  - r-gert=1.9.2=r41hf3f2ec2_0
  - r-getoptlong=1.0.5=r41hc72bb7e_1
  - r-ggplot2=3.4.2=r41hc72bb7e_0
  - r-ggrepel=0.9.3=r41h38f115c_0
  - r-ggridges=0.5.4=r41hc72bb7e_1
  - r-gh=1.4.0=r41hc72bb7e_0
  - r-gitcreds=0.1.2=r41hc72bb7e_1
  - r-globaloptions=0.1.2=r41ha770c72_1
  - r-globals=0.16.2=r41hc72bb7e_0
  - r-glue=1.6.2=r41h06615bd_1
  - r-goftest=1.2_3=r41h06615bd_1
  - r-googledrive=2.1.0=r41hc72bb7e_0
  - r-googlesheets4=1.1.0=r41h785f33e_0
  - r-gplots=3.1.3=r41hc72bb7e_1
  - r-gridextra=2.3=r41hc72bb7e_1004
  - r-gtable=0.3.3=r41hc72bb7e_0
  - r-gtools=3.9.4=r41h06615bd_0
  - r-harmony=1.1.0=r41h08d816e_0
  - r-haven=2.5.2=r41h38f115c_0
  - r-here=1.0.1=r41hc72bb7e_1
  - r-hexbin=1.28.3=r41hac0b197_0
  - r-highr=0.10=r41hc72bb7e_0
  - r-hms=1.1.3=r41hc72bb7e_0
  - r-htmltools=0.5.5=r41h38f115c_0
  - r-htmlwidgets=1.6.2=r41hc72bb7e_0
  - r-httpuv=1.6.11=r41ha503ecb_0
  - r-httr=1.4.6=r41hc72bb7e_0
  - r-httr2=0.2.3=r41hc72bb7e_0
  - r-ica=1.0_3=r41hc72bb7e_1
  - r-ids=1.0.1=r41hc72bb7e_2
  - r-igraph=1.4.2=r41h65ed38e_0
  - r-ini=0.3.1=r41hc72bb7e_1004
  - r-irdisplay=1.1=r41hd8ed1ab_1
  - r-irkernel=1.3.2=r41h785f33e_0
  - r-irlba=2.3.5.1=r41h5f7b363_0
  - r-isoband=0.2.7=r41h38f115c_1
  - r-iterators=1.0.14=r41hc72bb7e_1
  - r-jquerylib=0.1.4=r41hc72bb7e_1
  - r-jsonlite=1.8.5=r41h57805ef_0
  - r-kernsmooth=2.23_21=r41h13b3f57_0
  - r-knitr=1.43=r41hc72bb7e_0
  - r-labeling=0.4.2=r41hc72bb7e_2
  - r-lambda.r=1.2.4=r41hc72bb7e_2
  - r-later=1.3.1=r41ha503ecb_0
  - r-lattice=0.21_8=r41h133d619_0
  - r-lazyeval=0.2.2=r41h06615bd_3
  - r-leiden=0.4.3=r41hc72bb7e_1
  - r-lifecycle=1.0.3=r41hc72bb7e_1
  - r-listenv=0.9.0=r41hc72bb7e_0
  - r-lmtest=0.9_40=r41h8da6f51_1
  - r-lsei=1.3_0=r41hc3ea6d6_2
  - r-lubridate=1.9.2=r41h133d619_1
  - r-magrittr=2.0.3=r41h06615bd_1
  - r-mass=7.3_58.3=r41h133d619_0
  - r-matrix=1.5_4.1=r41h316c678_0
  - r-matrixstats=1.0.0=r41h57805ef_0
  - r-memoise=2.0.1=r41hc72bb7e_1
  - r-mgcv=1.8_42=r41he1ae0d6_0
  - r-mime=0.12=r41h06615bd_1
  - r-miniui=0.1.1.1=r41hc72bb7e_1003
  - r-modelr=0.1.11=r41hc72bb7e_0
  - r-munsell=0.5.0=r41hc72bb7e_1005
  - r-nabor=0.5.0=r41h7525677_1
  - r-nlme=3.1_162=r41hac0b197_0
  - r-npsurv=0.5_0=r41hc72bb7e_1
  - r-openssl=2.0.5=r41hb1dc35e_0
  - r-parallelly=1.36.0=r41hc72bb7e_0
  - r-patchwork=1.1.2=r41hc72bb7e_1
  - r-pbapply=1.7_0=r41hc72bb7e_0
  - r-pbdzmq=0.3_9=r41hfae1697_0
  - r-pillar=1.9.0=r41hc72bb7e_0
  - r-pkgbuild=1.4.0=r41hc72bb7e_0
  - r-pkgconfig=2.0.3=r41hc72bb7e_2
  - r-pkgdown=2.0.7=r41hc72bb7e_0
  - r-pkgload=1.3.2=r41hc72bb7e_0
  - r-plogr=0.2.0=r41hc72bb7e_1004
  - r-plotly=4.10.2=r41hc72bb7e_0
  - r-plyr=1.8.8=r41h7525677_0
  - r-png=0.1_8=r41h10cf519_0
  - r-polyclip=1.10_4=r41h7525677_0
  - r-powerlaw=0.70.6=r41hc72bb7e_2
  - r-pracma=2.4.2=r41hc72bb7e_1
  - r-praise=1.0.0=r41hc72bb7e_1006
  - r-presto=1.0.0=r41hdbdd923_1
  - r-prettyunits=1.1.1=r41hc72bb7e_2
  - r-processx=3.8.1=r41h133d619_0
  - r-profvis=0.3.8=r41h57805ef_0
  - r-progress=1.2.2=r41hc72bb7e_3
  - r-progressr=0.13.0=r41hc72bb7e_0
  - r-promises=1.2.0.1=r41h7525677_1
  - r-ps=1.7.5=r41h133d619_0
  - r-purrr=1.0.1=r41h133d619_0
  - r-r.methodss3=1.8.2=r41hc72bb7e_1
  - r-r.oo=1.25.0=r41hc72bb7e_1
  - r-r.utils=2.12.2=r41hc72bb7e_0
  - r-r6=2.5.1=r41hc72bb7e_1
  - r-ragg=1.2.5=r41hd65d3ba_0
  - r-rann=2.6.1=r41h7525677_3
  - r-rappdirs=0.3.3=r41h06615bd_1
  - r-rcmdcheck=1.4.0=r41h785f33e_1
  - r-rcolorbrewer=1.1_3=r41h785f33e_1
  - r-rcpp=1.0.10=r41h38f115c_0
  - r-rcppannoy=0.0.20=r41h7525677_0
  - r-rcpparmadillo=0.12.4.0.0=r41h08d816e_0
  - r-rcppeigen=0.3.3.9.3=r41h9f5de39_0
  - r-rcpphnsw=0.4.1=r41h7525677_1
  - r-rcppparallel=5.1.6=r41h38f115c_0
  - r-rcppprogress=0.4.2=r41hc72bb7e_2
  - r-rcpptoml=0.2.2=r41h38f115c_0
  - r-rcurl=1.98_1.12=r41h133d619_0
  - r-readr=2.1.4=r41h38f115c_0
  - r-readxl=1.4.2=r41h81ef4d7_0
  - r-rematch=1.0.1=r41hc72bb7e_1005
  - r-rematch2=2.1.2=r41hc72bb7e_2
  - r-remotes=2.4.2=r41hc72bb7e_1
  - r-repr=1.1.6=r41h785f33e_0
  - r-reprex=2.0.2=r41hc72bb7e_1
  - r-reshape2=1.4.4=r41h7525677_2
  - r-restfulr=0.0.15=r41h73dbb54_0
  - r-reticulate=1.30=r41ha503ecb_1
  - r-rgeos=0.6_3=r41h214e639_0
  - r-rhpcblasctl=0.23_42=r41h57805ef_0
  - r-rjson=0.2.21=r41h7525677_2
  - r-rlang=1.1.1=r41ha503ecb_0
  - r-rmarkdown=2.22=r41hc72bb7e_0
  - r-rocr=1.0_11=r41hc72bb7e_2
  - r-roxygen2=7.2.3=r41h38f115c_0
  - r-rprojroot=2.0.3=r41hc72bb7e_1
  - r-rspectra=0.16_1=r41h9f5de39_1
  - r-rsqlite=2.3.1=r41h38f115c_0
  - r-rstudioapi=0.14=r41hc72bb7e_1
  - r-rtsne=0.16=r41h37cf8d7_1
  - r-rversions=2.1.2=r41hc72bb7e_1
  - r-rvest=1.0.3=r41hc72bb7e_1
  - r-sass=0.4.6=r41ha503ecb_0
  - r-scales=1.2.1=r41hc72bb7e_1
  - r-scattermore=1.1=r41ha503ecb_0
  - r-sctransform=0.3.5=r41h9f5de39_1
  - r-selectr=0.4_2=r41hc72bb7e_2
  - r-sessioninfo=1.2.2=r41hc72bb7e_1
  - r-seurat=4.3.0=r41h38f115c_0
  - r-seuratobject=4.1.3=r41h38f115c_0
  - r-shape=1.4.6=r41ha770c72_1
  - r-shiny=1.7.4=r41h785f33e_0
  - r-sitmo=2.0.2=r41h7525677_1
  - r-snow=0.4_4=r41hc72bb7e_1
  - r-sourcetools=0.1.7_1=r41h38f115c_0
  - r-sp=1.6_1=r41h57805ef_0
  - r-spam=2.9_1=r41hb20cf53_1
  - r-spatstat.data=3.0_1=r41hc72bb7e_0
  - r-spatstat.explore=3.2_1=r41h57805ef_0
  - r-spatstat.geom=3.2_1=r41h57805ef_0
  - r-spatstat.random=3.1_5=r41ha503ecb_0
  - r-spatstat.sparse=3.0_1=r41h133d619_0
  - r-spatstat.univar=3.0_1=r41h2b5f3a1_1
  - r-spatstat.utils=3.1_0=r41h2b5f3a1_1
  - r-stringi=1.7.12=r41h1ae9187_0
  - r-stringr=1.5.0=r41h785f33e_0
  - r-survival=3.5_5=r41h133d619_0
  - r-sys=3.4.2=r41h57805ef_0
  - r-systemfonts=1.0.4=r41h0ff29ef_1
  - r-tensor=1.5=r41hc72bb7e_1004
  - r-testthat=3.1.8=r41ha503ecb_0
  - r-textshaping=0.3.6=r41hbb20487_4
  - r-tfmpvalue=0.0.9=r41h7525677_0
  - r-tibble=3.2.1=r41h133d619_1
  - r-tidyr=1.3.0=r41h38f115c_0
  - r-tidyselect=1.2.0=r41hc72bb7e_0
  - r-tidyverse=2.0.0=r41h785f33e_0
  - r-timechange=0.2.0=r41h38f115c_0
  - r-tinytex=0.45=r41hc72bb7e_0
  - r-tzdb=0.4.0=r41ha503ecb_0
  - r-urlchecker=1.0.1=r41hc72bb7e_1
  - r-usethis=2.2.0=r41hc72bb7e_0
  - r-utf8=1.2.3=r41h133d619_0
  - r-uuid=1.1_0=r41h06615bd_1
  - r-uwot=0.1.14=r41h7525677_1
  - r-vctrs=0.6.2=r41ha503ecb_0
  - r-viridislite=0.4.1=r41hc72bb7e_1
  - r-vroom=1.6.3=r41ha503ecb_0
  - r-waldo=0.5.1=r41hc72bb7e_0
  - r-whisker=0.4.1=r41hc72bb7e_0
  - r-withr=2.5.0=r41hc72bb7e_1
  - r-xfun=0.39=r41ha503ecb_0
  - r-xml=3.99_0.14=r41hb43fdd4_0
  - r-xml2=1.3.3=r41h044e5c7_2
  - r-xopen=1.0.0=r41hc72bb7e_1004
  - r-xtable=1.8_4=r41hc72bb7e_4
  - r-yaml=2.3.7=r41h133d619_0
  - r-zip=2.3.0=r41h133d619_0
  - r-zoo=1.8_12=r41h133d619_0
  - readline=8.2=h8228510_1
  - sed=4.8=he412f7d_0
  - sysroot_linux-64=2.17=h4a8ded7_17
  - tk=8.6.13=noxft_h4845f30_101
  - tktable=2.10=h8bc8fbc_6
  - tzdata=2024b=hc8b5060_0
  - xorg-kbproto=1.0.7=hb9d3cd8_1003
  - xorg-libice=1.0.10=h7f98852_0
  - xorg-libsm=1.2.3=hd9c2040_1000
  - xorg-libx11=1.8.4=h0b41bf4_0
  - xorg-libxau=1.0.11=hb9d3cd8_1
  - xorg-libxdmcp=1.1.5=hb9d3cd8_0
  - xorg-libxext=1.3.4=h0b41bf4_2
  - xorg-libxrender=0.9.10=h7f98852_1003
  - xorg-libxt=1.3.0=hd590300_0
  - xorg-renderproto=0.11.1=hb9d3cd8_1003
  - xorg-xextproto=7.3.0=hb9d3cd8_1004
  - xorg-xproto=7.0.31=hb9d3cd8_1008
  - xz=5.2.6=h166bdaf_0
  - zeromq=4.3.5=h59595ed_1
  - zlib=1.2.13=h4ab18f5_6
  - zstd=1.5.6=ha6fb4c9_0
prefix: /work/aaa/miniforge3/envs/archr-dev-r41
```
  
</details>

### SCAVENGE Pipeline
Software versions are documented in the [**cvejic-group/scavenge-smk-devmult**](https://github.com/cvejic-group/scavenge-smk-devmult) repository.

### Analysis
**`sessionInfo()` output**

<details>

```r
## R version 4.3.3 (2024-02-29)
## Platform: x86_64-conda-linux-gnu (64-bit)
## Running under: Ubuntu 22.04.4 LTS
## 
## Matrix products: default
## BLAS/LAPACK: /work/aaa/miniforge3/envs/scavenge_v102/lib/libopenblasp-r0.3.27.so;  LAPACK version 3.12.0
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## time zone: Europe/Berlin
## tzcode source: system (glibc)
## 
## attached base packages:
## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] pheatmap_1.0.12             magrittr_2.0.3             
##  [3] lubridate_1.9.3             forcats_1.0.0              
##  [5] stringr_1.5.1               dplyr_1.1.4                
##  [7] purrr_1.0.2                 readr_2.1.5                
##  [9] tidyr_1.3.1                 tibble_3.2.1               
## [11] ggplot2_3.5.1               tidyverse_2.0.0            
## [13] HDF5Array_1.30.0            rhdf5_2.46.1               
## [15] DelayedArray_0.28.0         SparseArray_1.2.2          
## [17] S4Arrays_1.2.0              abind_1.4-5                
## [19] Matrix_1.6-5                SummarizedExperiment_1.32.0
## [21] Biobase_2.62.0              GenomicRanges_1.54.1       
## [23] GenomeInfoDb_1.38.1         IRanges_2.36.0             
## [25] S4Vectors_0.40.2            BiocGenerics_0.48.1        
## [27] MatrixGenerics_1.14.0       matrixStats_1.4.1          
## 
## loaded via a namespace (and not attached):
##  [1] gtable_0.3.5            xfun_0.48               bslib_0.8.0            
##  [4] lattice_0.22-6          tzdb_0.4.0              rhdf5filters_1.14.1    
##  [7] vctrs_0.6.5             tools_4.3.3             bitops_1.0-9           
## [10] generics_0.1.3          parallel_4.3.3          fansi_1.0.6            
## [13] highr_0.11              pkgconfig_2.0.3         RColorBrewer_1.1-3     
## [16] lifecycle_1.0.4         GenomeInfoDbData_1.2.11 farver_2.1.2           
## [19] compiler_4.3.3          munsell_0.5.1           htmltools_0.5.8.1      
## [22] sass_0.4.9              RCurl_1.98-1.16         yaml_2.3.10            
## [25] pillar_1.9.0            crayon_1.5.3            jquerylib_0.1.4        
## [28] cachem_1.1.0            tidyselect_1.2.1        digest_0.6.37          
## [31] stringi_1.8.4           fastmap_1.2.0           grid_4.3.3             
## [34] colorspace_2.1-1        cli_3.6.3               utf8_1.2.4             
## [37] withr_3.0.1             scales_1.3.0            bit64_4.5.2            
## [40] timechange_0.3.0        rmarkdown_2.28          XVector_0.42.0         
## [43] bit_4.5.0               hms_1.1.3               evaluate_1.0.0         
## [46] knitr_1.48              viridisLite_0.4.2       rlang_1.1.4            
## [49] glue_1.8.0              vroom_1.6.5             rstudioapi_0.16.0      
## [52] jsonlite_1.8.9          R6_2.5.1                Rhdf5lib_1.24.0        
## [55] zlibbioc_1.48.0
```
  
</details>

**Conda Environment**

<details>

```yaml
name: scavenge_v102
channels:
  - merv
  - bioconda
  - conda-forge
dependencies:
  - _libgcc_mutex=0.1=conda_forge
  - _openmp_mutex=4.5=2_gnu
  - _r-mutex=1.0.1=anacondar_1
  - argcomplete=3.5.1=pyhd8ed1ab_0
  - binutils_impl_linux-64=2.43=h4bf12b8_1
  - bioconductor-annotate=1.80.0=r43hdfd78af_0
  - bioconductor-annotationdbi=1.64.1=r43hdfd78af_0
  - bioconductor-biobase=2.62.0=r43ha9d7317_2
  - bioconductor-biocgenerics=0.48.1=r43hdfd78af_2
  - bioconductor-biocio=1.12.0=r43hdfd78af_0
  - bioconductor-biocparallel=1.36.0=r43hf17093f_2
  - bioconductor-biostrings=2.70.1=r43ha9d7317_2
  - bioconductor-bsgenome=1.70.1=r43hdfd78af_0
  - bioconductor-bsgenome.hsapiens.ucsc.hg19=1.4.3=r43hdfd78af_8
  - bioconductor-bsgenome.hsapiens.ucsc.hg38=1.4.5=r43hdfd78af_2
  - bioconductor-chromvar=1.24.0=r43hf17093f_0
  - bioconductor-cner=1.38.0=r43ha9d7317_1
  - bioconductor-data-packages=20231203=hdfd78af_0
  - bioconductor-delayedarray=0.28.0=r43ha9d7317_2
  - bioconductor-dirichletmultinomial=1.44.0=r43hee7dd41_1
  - bioconductor-genomeinfodb=1.38.1=r43hdfd78af_1
  - bioconductor-genomeinfodbdata=1.2.11=r43hdfd78af_1
  - bioconductor-genomicalignments=1.38.0=r43ha9d7317_1
  - bioconductor-genomicranges=1.54.1=r43ha9d7317_2
  - bioconductor-go.db=3.18.0=r43hdfd78af_0
  - bioconductor-hdf5array=1.30.0=r43ha9d7317_1
  - bioconductor-iranges=2.36.0=r43ha9d7317_2
  - bioconductor-keggrest=1.42.0=r43hdfd78af_0
  - bioconductor-matrixgenerics=1.14.0=r43hdfd78af_3
  - bioconductor-rhdf5=2.46.1=r43hf17093f_1
  - bioconductor-rhdf5filters=1.14.1=r43hf17093f_1
  - bioconductor-rhdf5lib=1.24.0=r43ha9d7317_2
  - bioconductor-rhtslib=2.4.0=r43ha9d7317_2
  - bioconductor-rsamtools=2.18.0=r43hf17093f_2
  - bioconductor-rtracklayer=1.62.0=r43ha9d7317_1
  - bioconductor-s4arrays=1.2.0=r43ha9d7317_2
  - bioconductor-s4vectors=0.40.2=r43ha9d7317_2
  - bioconductor-seqlogo=1.68.0=r43hdfd78af_0
  - bioconductor-sparsearray=1.2.2=r43ha9d7317_2
  - bioconductor-summarizedexperiment=1.32.0=r43hdfd78af_0
  - bioconductor-tfbstools=1.40.0=r43ha9d7317_1
  - bioconductor-xvector=0.42.0=r43ha9d7317_2
  - bioconductor-zlibbioc=1.48.0=r43ha9d7317_2
  - bwidget=1.9.14=ha770c72_1
  - bzip2=1.0.8=h4bc722e_7
  - c-ares=1.33.1=heb4867d_0
  - ca-certificates=2024.8.30=hbcca054_0
  - cairo=1.18.0=hebfffa5_3
  - curl=8.10.1=hbbe4b11_0
  - expat=2.6.3=h5888daf_0
  - font-ttf-dejavu-sans-mono=2.37=hab24e00_0
  - font-ttf-inconsolata=3.000=h77eed37_0
  - font-ttf-source-code-pro=2.038=h77eed37_0
  - font-ttf-ubuntu=0.83=h77eed37_3
  - fontconfig=2.14.2=h14ed4e7_0
  - fonts-conda-ecosystem=1=0
  - fonts-conda-forge=1=0
  - freetype=2.12.1=h267a509_2
  - fribidi=1.0.10=h36c2ea0_0
  - gcc_impl_linux-64=14.1.0=h3c94d91_1
  - gfortran_impl_linux-64=14.1.0=he4a1faa_1
  - glpk=5.0=h445213a_0
  - gmp=6.3.0=hac33072_2
  - graphite2=1.3.13=h59595ed_1003
  - gsl=2.7=he838d99_0
  - gxx_impl_linux-64=14.1.0=h8d00ecb_1
  - harfbuzz=9.0.0=hda332d3_1
  - icu=75.1=he02047a_0
  - jq=1.7.1=hd590300_0
  - kernel-headers_linux-64=3.10.0=he073ed8_17
  - keyutils=1.6.1=h166bdaf_0
  - krb5=1.21.3=h659f571_0
  - ld_impl_linux-64=2.43=h712a8e2_1
  - lerc=4.0.0=h27087fc_0
  - libblas=3.9.0=24_linux64_openblas
  - libcblas=3.9.0=24_linux64_openblas
  - libcurl=8.10.1=hbbe4b11_0
  - libdeflate=1.22=hb9d3cd8_0
  - libedit=3.1.20191231=he28a2e2_2
  - libev=4.33=hd590300_2
  - libexpat=2.6.3=h5888daf_0
  - libffi=3.4.2=h7f98852_5
  - libgcc=14.1.0=h77fa898_1
  - libgcc-devel_linux-64=14.1.0=h5d3d1c9_101
  - libgcc-ng=14.1.0=h69a702a_1
  - libgfortran=14.1.0=h69a702a_1
  - libgfortran-ng=14.1.0=h69a702a_1
  - libgfortran5=14.1.0=hc5f4f2c_1
  - libglib=2.82.1=h2ff4ddf_0
  - libgomp=14.1.0=h77fa898_1
  - libiconv=1.17=hd590300_2
  - libjpeg-turbo=3.0.0=hd590300_1
  - liblapack=3.9.0=24_linux64_openblas
  - libmpdec=4.0.0=h4bc722e_0
  - libnghttp2=1.58.0=h47da74e_1
  - libopenblas=0.3.27=pthreads_hac2b453_1
  - libpng=1.6.44=hadc24fc_0
  - libsanitizer=14.1.0=hcba0ae0_1
  - libsqlite=3.46.1=hadc24fc_0
  - libssh2=1.11.0=h0841786_0
  - libstdcxx=14.1.0=hc0a3c3a_1
  - libstdcxx-devel_linux-64=14.1.0=h5d3d1c9_101
  - libstdcxx-ng=14.1.0=h4852527_1
  - libtiff=4.7.0=he137b08_1
  - libuuid=2.38.1=h0b41bf4_0
  - libuv=1.49.0=hb9d3cd8_0
  - libwebp-base=1.4.0=hd590300_0
  - libxcb=1.17.0=h8a09558_0
  - libxml2=2.12.7=he7c6b58_4
  - libzlib=1.3.1=hb9d3cd8_2
  - make=4.4.1=hb9d3cd8_2
  - ncurses=6.5=he02047a_1
  - oniguruma=6.9.9=hd590300_0
  - openssl=3.3.2=hb9d3cd8_0
  - pandoc=3.5=ha770c72_0
  - pango=1.54.0=h4c5309f_1
  - pcre2=10.44=hba22ea6_2
  - pip=24.2=pyh145f28c_1
  - pixman=0.43.2=h59595ed_0
  - pthread-stubs=0.4=hb9d3cd8_1002
  - python=3.13.0=h9ebbce0_100_cp313
  - python_abi=3.13=5_cp313
  - pyyaml=6.0.2=py313h536fd9c_1
  - r-abind=1.4_5=r43hc72bb7e_1006
  - r-askpass=1.2.1=r43h2b5f3a1_0
  - r-assertthat=0.2.1=r43hc72bb7e_5
  - r-backports=1.5.0=r43hb1dbf0f_1
  - r-base=4.3.3=h68df6d5_15
  - r-base64enc=0.1_3=r43hb1dbf0f_1007
  - r-beeswarm=0.4.0=r43hdb488b9_4
  - r-bh=1.84.0_0=r43hc72bb7e_1
  - r-bit=4.5.0=r43h2b5f3a1_0
  - r-bit64=4.5.2=r43h2b5f3a1_0
  - r-bitops=1.0_9=r43h2b5f3a1_0
  - r-blob=1.2.4=r43hc72bb7e_2
  - r-broom=1.0.7=r43hc72bb7e_0
  - r-bslib=0.8.0=r43hc72bb7e_0
  - r-buencolors=0.5.6=r43hdfd78af_0
  - r-cachem=1.1.0=r43hb1dbf0f_1
  - r-cairo=1.6_2=r43h232ff4d_2
  - r-callr=3.7.6=r43hc72bb7e_1
  - r-catools=1.18.3=r43h93ab643_0
  - r-cellranger=1.1.0=r43hc72bb7e_1007
  - r-cli=3.6.3=r43h0d4f4ea_1
  - r-clipr=0.8.0=r43hc72bb7e_3
  - r-codetools=0.2_20=r43hc72bb7e_1
  - r-colorspace=2.1_1=r43hdb488b9_0
  - r-commonmark=1.9.2=r43h2b5f3a1_0
  - r-conflicted=1.2.0=r43h785f33e_2
  - r-cpp11=0.5.0=r43hc72bb7e_0
  - r-crayon=1.5.3=r43hc72bb7e_1
  - r-crosstalk=1.2.1=r43hc72bb7e_1
  - r-curl=5.2.1=r43h6b349a7_1
  - r-data.table=1.15.4=r43h5f06984_1
  - r-dbi=1.2.3=r43hc72bb7e_1
  - r-dbplyr=2.5.0=r43hc72bb7e_1
  - r-digest=0.6.37=r43h0d4f4ea_0
  - r-dplyr=1.1.4=r43h0d4f4ea_1
  - r-dt=0.33=r43hc72bb7e_1
  - r-dtplyr=1.3.1=r43hc72bb7e_2
  - r-ellipsis=0.3.2=r43hb1dbf0f_3
  - r-evaluate=1.0.0=r43hc72bb7e_0
  - r-fansi=1.0.6=r43hb1dbf0f_1
  - r-farver=2.1.2=r43ha18555a_1
  - r-fastmap=1.2.0=r43ha18555a_1
  - r-fontawesome=0.5.2=r43hc72bb7e_1
  - r-forcats=1.0.0=r43hc72bb7e_2
  - r-formatr=1.14=r43hc72bb7e_2
  - r-fs=1.6.4=r43ha18555a_1
  - r-futile.logger=1.4.3=r43hc72bb7e_1006
  - r-futile.options=1.0.1=r43hc72bb7e_1005
  - r-gargle=1.5.2=r43h785f33e_1
  - r-gchromvar=0.3.2=r43hdfd78af_0
  - r-generics=0.1.3=r43hc72bb7e_3
  - r-ggbeeswarm=0.7.2=r43hc72bb7e_2
  - r-ggplot2=3.5.1=r43hc72bb7e_1
  - r-ggrastr=1.0.2=r43hc72bb7e_2
  - r-glue=1.8.0=r43h2b5f3a1_0
  - r-googledrive=2.1.1=r43hc72bb7e_2
  - r-googlesheets4=1.1.1=r43h785f33e_2
  - r-gtable=0.3.5=r43hc72bb7e_1
  - r-gtools=3.9.5=r43hb1dbf0f_1
  - r-haven=2.5.4=r43h0d4f4ea_1
  - r-hexbin=1.28.4=r43hb67ce94_0
  - r-highr=0.11=r43hc72bb7e_1
  - r-hms=1.1.3=r43hc72bb7e_2
  - r-htmltools=0.5.8.1=r43ha18555a_1
  - r-htmlwidgets=1.6.4=r43h785f33e_3
  - r-httpuv=1.6.15=r43ha18555a_1
  - r-httr=1.4.7=r43hc72bb7e_1
  - r-ids=1.0.1=r43hc72bb7e_4
  - r-igraph=2.0.3=r43h510f1ce_1
  - r-irlba=2.3.5.1=r43h0d28552_3
  - r-isoband=0.2.7=r43ha18555a_3
  - r-jquerylib=0.1.4=r43hc72bb7e_3
  - r-jsonlite=1.8.9=r43h2b5f3a1_0
  - r-knitr=1.48=r43hc72bb7e_0
  - r-labeling=0.4.3=r43hc72bb7e_1
  - r-lambda.r=1.2.4=r43hc72bb7e_4
  - r-later=1.3.2=r43ha18555a_1
  - r-lattice=0.22_6=r43hb1dbf0f_1
  - r-lazyeval=0.2.2=r43hb1dbf0f_5
  - r-lifecycle=1.0.4=r43hc72bb7e_1
  - r-lubridate=1.9.3=r43hdb488b9_1
  - r-magrittr=2.0.3=r43hb1dbf0f_3
  - r-mass=7.3_60.0.1=r43hb1dbf0f_1
  - r-matrix=1.6_5=r43he966344_1
  - r-matrixstats=1.4.1=r43h2b5f3a1_0
  - r-memoise=2.0.1=r43hc72bb7e_3
  - r-mervdown=0.2.0=r43_0
  - r-mgcv=1.9_1=r43h0d28552_1
  - r-mime=0.12=r43hb1dbf0f_3
  - r-miniui=0.1.1.1=r43hc72bb7e_1005
  - r-modelr=0.1.11=r43hc72bb7e_2
  - r-munsell=0.5.1=r43hc72bb7e_1
  - r-nabor=0.5.0=r43h0d4f4ea_3
  - r-nlme=3.1_165=r43hbcb9c34_1
  - r-openssl=2.2.2=r43he8289e2_0
  - r-pheatmap=1.0.12=r43hc72bb7e_5
  - r-pillar=1.9.0=r43hc72bb7e_2
  - r-pkgconfig=2.0.3=r43hc72bb7e_4
  - r-plogr=0.2.0=r43hc72bb7e_1006
  - r-plotly=4.10.4=r43hc72bb7e_1
  - r-plyr=1.8.9=r43ha18555a_1
  - r-png=0.1_8=r43h21f035c_2
  - r-powerlaw=0.80.0=r43hc72bb7e_1
  - r-pracma=2.4.4=r43hc72bb7e_1
  - r-prettyunits=1.2.0=r43hc72bb7e_1
  - r-processx=3.8.4=r43hb1dbf0f_1
  - r-progress=1.2.3=r43hc72bb7e_1
  - r-promises=1.3.0=r43ha18555a_1
  - r-ps=1.7.7=r43hdb488b9_0
  - r-purrr=1.0.2=r43hdb488b9_1
  - r-r.methodss3=1.8.2=r43hc72bb7e_3
  - r-r.oo=1.26.0=r43hc72bb7e_1
  - r-r.utils=2.12.3=r43hc72bb7e_1
  - r-r6=2.5.1=r43hc72bb7e_3
  - r-ragg=1.3.3=r43h9aa3752_0
  - r-rann=2.6.2=r43h93ab643_0
  - r-rappdirs=0.3.3=r43hb1dbf0f_3
  - r-rcolorbrewer=1.1_3=r43h785f33e_3
  - r-rcpp=1.0.13=r43h0d4f4ea_0
  - r-rcpparmadillo=14.0.2_1=r43hb79369c_0
  - r-rcppeigen=0.3.4.0.2=r43hb79369c_0
  - r-rcurl=1.98_1.16=r43he8228da_1
  - r-readr=2.1.5=r43h0d4f4ea_1
  - r-readxl=1.4.3=r43he58e087_1
  - r-rematch=2.0.0=r43hc72bb7e_1
  - r-rematch2=2.1.2=r43hc72bb7e_4
  - r-remotes=2.5.0=r43hc72bb7e_1
  - r-reprex=2.1.1=r43hc72bb7e_1
  - r-reshape2=1.4.4=r43h0d4f4ea_4
  - r-restfulr=0.0.15=r43h56115f1_4
  - r-rjson=0.2.23=r43h93ab643_0
  - r-rlang=1.1.4=r43ha18555a_1
  - r-rmarkdown=2.28=r43hc72bb7e_0
  - r-rsqlite=2.3.7=r43h0d4f4ea_0
  - r-rstudioapi=0.16.0=r43hc72bb7e_1
  - r-rtsne=0.17=r43h4387864_1
  - r-rvest=1.0.4=r43hc72bb7e_1
  - r-sass=0.4.9=r43ha18555a_1
  - r-scales=1.3.0=r43hc72bb7e_1
  - r-scavenge=1.0.2=r43hdfd78af_0
  - r-selectr=0.4_2=r43hc72bb7e_4
  - r-shiny=1.9.1=r43h785f33e_0
  - r-snow=0.4_4=r43hc72bb7e_3
  - r-sourcetools=0.1.7_1=r43ha18555a_2
  - r-stringi=1.8.4=r43h33cde33_3
  - r-stringr=1.5.1=r43h785f33e_1
  - r-sys=3.4.3=r43h2b5f3a1_0
  - r-systemfonts=1.1.0=r43h38d38ca_1
  - r-textshaping=0.4.0=r43ha47bcaa_2
  - r-tfmpvalue=0.0.9=r43h0d4f4ea_2
  - r-tibble=3.2.1=r43hdb488b9_3
  - r-tidyr=1.3.1=r43h0d4f4ea_1
  - r-tidyselect=1.2.1=r43hc72bb7e_1
  - r-tidyverse=2.0.0=r43h785f33e_2
  - r-timechange=0.3.0=r43ha18555a_1
  - r-tinytex=0.53=r43hc72bb7e_0
  - r-tzdb=0.4.0=r43ha18555a_2
  - r-utf8=1.2.4=r43hb1dbf0f_1
  - r-uuid=1.2_1=r43hdb488b9_0
  - r-vctrs=0.6.5=r43h0d4f4ea_1
  - r-vipor=0.4.7=r43hc72bb7e_1
  - r-viridislite=0.4.2=r43hc72bb7e_2
  - r-vroom=1.6.5=r43h0d4f4ea_1
  - r-withr=3.0.1=r43hc72bb7e_0
  - r-xfun=0.48=r43h93ab643_0
  - r-xml=3.99_0.17=r43he716329_1
  - r-xml2=1.3.6=r43h8194278_2
  - r-xtable=1.8_4=r43hc72bb7e_6
  - r-yaml=2.3.10=r43hdb488b9_0
  - readline=8.2=h8228510_1
  - sed=4.8=he412f7d_0
  - setuptools=75.1.0=pyhd8ed1ab_0
  - sysroot_linux-64=2.17=h4a8ded7_17
  - tk=8.6.13=noxft_h4845f30_101
  - tktable=2.10=h8bc8fbc_6
  - toml=0.10.2=pyhd8ed1ab_0
  - tomlkit=0.13.2=pyha770c72_0
  - tzdata=2024b=hc8b5060_0
  - xmltodict=0.14.0=pyhd8ed1ab_0
  - xorg-libice=1.1.1=hb9d3cd8_1
  - xorg-libsm=1.2.4=he73a12e_1
  - xorg-libx11=1.8.10=h4f16b4b_0
  - xorg-libxau=1.0.11=hb9d3cd8_1
  - xorg-libxdmcp=1.1.5=hb9d3cd8_0
  - xorg-libxext=1.3.6=hb9d3cd8_0
  - xorg-libxrender=0.9.11=hb9d3cd8_1
  - xorg-libxt=1.3.0=hb9d3cd8_2
  - xorg-xorgproto=2024.1=hb9d3cd8_1
  - xz=5.2.6=h166bdaf_0
  - yaml=0.2.5=h7f98852_2
  - yq=3.4.3=pyhd8ed1ab_0
  - zlib=1.3.1=hb9d3cd8_2
  - zstd=1.5.6=ha6fb4c9_0
prefix: /work/aaa/miniforge3/envs/scavenge_v102
```

</details>
