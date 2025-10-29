#!/usr/bin/env bash

# WORK DIR
WORKDIR=$HOME/project/20231127_DevM/sinto/bam_per_lib
PPN=30
libs=$1

# Tools
tool_dir=$HOME/software

# conda
source $tool_dir/anaconda3/bin/activate
conda activate sinto

for lib in "${libs[@]}"
do
  bam=/work/home/data/cellranger-arc/cellranger-arc202_count_${lib}_GRCh38-2020-A-2_0_0/atac_possorted_bam.bam
  #sinto filterbarcodes -p $PPN -b $bam -c $WORKDIR/$lib/cells.tsv --outdir $WORKDIR/$lib
  sinto filterbarcodes -p $PPN -b $bam -c $WORKDIR/$lib/cells_hsc.tsv --outdir $WORKDIR/$lib
done


