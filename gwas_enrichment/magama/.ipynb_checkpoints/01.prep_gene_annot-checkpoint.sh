#!/usr/bin/env bash

# usage
# bash 01.prep_gene_annot.sh celltype

work_dir=~/project/20231127_DevM/gwas/MAGMA_gene
out_dir=$work_dir/gene

# cell from input
cell=$1

# magma
tool_dir=$HOME/software
path2magma=$tool_dir/magma_v1.10/magma
magma_files=/work/DevM_analysis/data/magma_files

# peak-gene links
pg_bedpe=/work/DevM_analysis/data/E2G/lmmPG/$cell.lmmPG_FDR20.bedpe

# to bed
cut -f 1-3,7 $pg_bedpe | sed 's/:/\t/g' | cut -f 1-3,5 > tmp.bed

# liftover to hg19
liftOver tmp.bed ~/RefData/ucsc/liftOver/hg38ToHg19.over.chain.gz tmp2.bed tmp.unlifted

# intersect with SNP
sed 's/^chr//g' tmp2.bed | sort -k1,1V -k2,2n > $out_dir/${cell}.cCRE_hg19.bed
bedtools intersect -a $magma_files/g1000_eur/g1000_eur.bed4 -b $out_dir/${cell}.cCRE_hg19.bed -wa -wb > $out_dir/${cell}.snp_ccre_gene.txt

# Create a MAGMA gene annotation file
python make_genes_annot.py \
  --intersect $out_dir/${cell}.snp_ccre_gene.txt \
  --gene-pos ~/RefData/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes_pos.txt \
  --out $out_dir/${cell}.genes.annot


