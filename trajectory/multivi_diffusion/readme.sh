#!/usr/bin/env bash

conda activate scFates

# script to generate mvi space, B lineage as example
python 01.prep_atac.py
Rscript 01.prep_rna.R
python 02.prep_wnn.py
python 03.multivi.py

# diffusion map
# blood
python diffusion_mvi.py /work/DevM_analysis/01.annotation/11.subclustering/blood/data/FL_multiVI.h5ad FL_blood_diffusion > logs/FL_blood_diffusion.log 2>&1

# B
for pcw in `seq 5 18`
do
  python diffusion_mvi.py /work/DevM_analysis/01.annotation/11.subclustering/BlineageByPCW/data/PCW${pcw}_multiVI.h5ad Blineage_diffusion_PCW${pcw} > logs/Blineage_diffusion_PCW${pcw}.log 2>&1
done

# MEMP
for pcw in `seq 5 18`
do
  python diffusion_mvi.py /work/DevM_analysis/01.annotation/11.subclustering/ErythroByPCW/data/PCW${pcw}_multiVI.h5ad Erythro_diffusion_PCW${pcw} > logs/Erythro_diffusion_PCW${pcw}.log 2>&1
done

# Myelo
for pcw in `seq 5 18`
do
  python diffusion_mvi.py /work/DevM_analysis/01.annotation/11.subclustering/MyeloByPCW/data/PCW${pcw}_multiVI.h5ad Myelo_diffusion_PCW${pcw} > logs/Myelo_diffusion_PCW${pcw}.log 2>&1
done
