#!/usr/bin/env bash

# usage:
# bash 02.run_magma.sh cell pheno N

work_dir=~/project/20231127_DevM/gwas
out_dir=$work_dir/MAGMA_gene/result

# cell from input
cell=$1
pheno=$2
N=$3

# magma
tool_dir=$HOME/software
path2magma=$tool_dir/magma_v1.10/magma
magma_files=/work/DevM_analysis/data/magma_files

# Run MAGMA gene analysis
$path2magma \
  --bfile $magma_files/g1000_eur/g1000_eur synonyms=$magma_files/dbsnp151.synonyms/dbsnp151.synonyms \
  --pval $work_dir/MAGMA_cCRE/snpp/${pheno}.tsv N=$N \
  --gene-annot $work_dir/MAGMA_gene/gene/${cell}.genes.annot \
  --out $out_dir/${pheno}_${cell}_cCRE_analysis

