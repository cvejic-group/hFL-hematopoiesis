#!/usr/bin/env bash

# Argument parsing
while getopts "i:b:k:o:t:" opt; do
  case "$opt" in
    i) bam=${OPTARG};;
    b) barcode=${OPTARG};;
    k) cluster=${OPTARG};;
    o) outdir=${OPTARG};;
    t) PPN=${OPTARG};;
    *) echo -e "\033[1;31mERROR: invalid option, only accept -i, -b, -k, -o or -t.\033[0m" 1>&2; exit 1;;
  esac
done

# Tools
tool_dir=$HOME/software
path2vartrix=$tool_dir/vartrix_linux

# Anno
anno_dir=$HOME/RefData
ref_genome=$anno_dir/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa
ref_vcf=$anno_dir/souporcell/common_variants_hg38.vcf

# conda
source $tool_dir/anaconda3/bin/activate
conda activate souporcell
export PATH=$HOME/software/souporcell:$PATH
export PATH=$HOME/software/souporcell/souporcell/target/release:$PATH
export PATH=$HOME/software/souporcell/troublet/target/release:$PATH # for doublets

# Run
if [ ! -d $outdir/$cluster ]; then mkdir -p $outdir/$cluster; fi
souporcell_pipeline.py -i $bam -b $barcode -f $ref_genome --common_variants $ref_vcf -t $PPN -o $outdir/$cluster -k $cluster --skip_remap SKIP_REMAP

# run k=1/2/3
if [ "$cluster"x != "1"x ]; then
  if [ ! -d $outdir/1 ]; then mkdir -p $outdir/1; fi
  souporcell_pipeline.py -i $bam -b $barcode -f $ref_genome --common_variants $ref_vcf -t $PPN -o $outdir/1 -k 1 --skip_remap SKIP_REMAP
fi
if [ "$cluster"x != "2"x ]; then
  if [ ! -d $outdir/2 ]; then mkdir -p $outdir/2; fi
  souporcell_pipeline.py -i $bam -b $barcode -f $ref_genome --common_variants $ref_vcf -t $PPN -o $outdir/2 -k 2 --skip_remap SKIP_REMAP
fi
if [ "$cluster"x != "3"x ]; then
  if [ ! -d $outdir/3 ]; then mkdir -p $outdir/3; fi
  souporcell_pipeline.py -i $bam -b $barcode -f $ref_genome --common_variants $ref_vcf -t $PPN -o $outdir/3 -k 3 --skip_remap SKIP_REMAP
fi





