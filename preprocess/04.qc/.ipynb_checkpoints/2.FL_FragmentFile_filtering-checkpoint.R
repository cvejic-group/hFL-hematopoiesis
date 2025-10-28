# run
# cd /work/DevM_analysis/01.annotation/02.cleandata/
# nohup Rscript 02.2.FL_FragmentFile_filtering.R > 02.2.FL.log & 

.libPaths("/work/home/Software/R/")

library(Seurat) 
library(Signac)
library(dplyr)

meta_cells <- read.table("/work/DevM_analysis/01.annotation/02.cleandata/data/FL_cell-info_filtering.txt", sep = "\t")
sam_info <- read.table("/work/DevM_analysis/utils/sam_info.DevM.12.08.24.txt", sep = "\t", header = T)

##---------------------------------------------##
##---------3. Fragment file filtering----------##
##---------------------------------------------##

for(s in unique(sam_info$libraryID)){
  fpath <- paste0("/work/cellranger-arc/", 
                  (sam_info %>% filter(libraryID == s))[1, "CellRanger"], 
                  "/atac_fragments.tsv.gz")
  outpath <- paste0("/work/DevM_analysis/01.annotation/02.cleandata/data/filtered_frag_files/", s, "_HighQualityCells_atac_fragments.tsv.gz")
  
  cat(s, "\n")
  FilterCells(
    fragments = fpath,
    cells = (meta_cells %>% filter(HighQualityCell == 1) %>% filter(libraryID == s))[,"gex_barcode"],
    outfile = outpath
  )
}
