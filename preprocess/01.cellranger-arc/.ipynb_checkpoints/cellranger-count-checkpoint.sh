#!/usr/bin/env bash

# WORK DIR
WORKDIR=$HOME/Data_processing

# PPN
PPN=16

# dataset and sample
dataset=$1
sample=$2

# Anno
ANNODIR=$HOME/RefData

# tools dir
TOOLDIR=/work/home/software

# out (cellranger wants to create the folder)
output=$WORKDIR/cellranger/$dataset/$sample

# force to analyze the GEX part of ARC
$TOOLDIR/cellranger-8.0.1/cellranger count --localcores=$PPN --create-bam=true --id=$sample --fastqs=$WORKDIR/rawfastq/$dataset/${sample} --transcriptome=$ANNODIR/refdata-gex-GRCh38-2024-A --output-dir=$output --chemistry=ARC-v1


