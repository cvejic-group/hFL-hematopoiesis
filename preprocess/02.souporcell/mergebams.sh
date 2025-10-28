#!/usr/bin/env bash

work_dir=/work/home/project/20231127_DevM/souporcell/mergebams
PPN=24

# cr dir
cr_dir=/work/home/data/cellranger-arc

# Tools
tool_dir=/work/home/software

# Anno
anno_dir=/work/home/RefData

# lib1 and lib2
donor=$1
sam1=lib1; sam2=lib2
cr1=cellranger-arc202_count_${sam1}_GRCh38-2020-A-2_0_0
cr2=cellranger-arc202_count_${sam2}_GRCh38-2020-A-2_0_0

# run
bams=$cr_dir/$cr1/gex_possorted_bam.bam,$cr_dir/$cr2/gex_possorted_bam.bam
barcodes=$cr_dir/$cr1/filtered_feature_bc_matrix/barcodes.tsv.gz,$cr_dir/$cr2/filtered_feature_bc_matrix/barcodes.tsv.gz
labels=${sam1}_,${sam2}_
outdir=$work_dir/$donor
if [ ! -d $outdir ]; then mkdir -p $outdir; fi
$tool_dir/mergebams_linux -i $bams -b $barcodes -l $labels -t $PPN -o $outdir

# sort bam
$tool_dir/samtools/1.18/samtools sort -@ $PPN -o $outdir/out_bam.sorted.bam $outdir/out_bam.bam
# build index
$tool_dir/samtools/1.18/samtools index -@ $PPN $outdir/out_bam.sorted.bam



