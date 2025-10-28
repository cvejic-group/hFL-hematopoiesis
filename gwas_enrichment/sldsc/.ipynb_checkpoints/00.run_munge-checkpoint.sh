#!/usr/bin/env bash

# nohup bash 00.run_munge.sh >> 00.run_munge.log &

work_dir=~/project/20231127_DevM/gwas
data_dir=/work/DevM_analysis/data/gwas_results/MungeSumstats_GRCh37/ldsc_munge_ready
out_dir=$work_dir/sLDSC/munge_sumstats

# Tools
tool_dir=$HOME/software
path2ldsc=$tool_dir/ldsc

# Anno
ldsc_files=/work/DevM_analysis/data/S-LDSC_reference_files

# conda
source $tool_dir/anaconda3/bin/activate
conda activate ldsc

# traits
df_in=$work_dir/gwas_traits.tsv

# run
sed -n '2,$p' $df_in | cut -f 2 | while read pheno
do
  echo $pheno
  python $path2ldsc/munge_sumstats.py \
    --chunksize 500000 \
    --sumstats $data_dir/$pheno.tsv.gz \
    --merge-alleles $ldsc_files/w_hm3.snplist \
    --out $out_dir/$pheno
done

