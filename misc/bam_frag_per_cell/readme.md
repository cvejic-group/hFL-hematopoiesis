
## About

prep BAMs and fragment files per cell and per cell per PCW


## How

* `01.prep_cells.R`: prep barcodes of each type

* `02.filterbarcodes.sh`: use sinto to separate ATAC BAM file of each library into different cell types (output `bam_per_lib`)

* `03.prep_mergebam.py` and `03.merge_libs.sh`: merge BAM files of each library by cell type and by PCW  (output `bam_per_cell`)

* `04.prep_dedup.py` and `04.dedup.sh`: run de-dup (output `bam_per_cell`)

* `05.prep_fragments.py` and `05.prep_frag.sh`: generate fragment files of cell type and each PCW (output `frag_per_cell`)



