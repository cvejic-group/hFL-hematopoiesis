#!/usr/bin/env bash

# nohup bash 03.peak_ldscore.sh > 03.peak_ldscore.log &

# constrain threads usage by ldsc
export OMP_NUM_THREADS=10
export OPENBLAS_NUM_THREADS=10
export MKL_NUM_THREADS=10
export NUMEXPR_NUM_THREADS=10

work_dir=~/project/20231127_DevM/gwas/sLDSC
peak_dir=$work_dir/atac_peaks
d=FDR5_LFC1_padding_1k

# ldsc
tool_dir=$HOME/software
path2ldsc=$tool_dir/ldsc
ldsc_files=/work/DevM_analysis/data/S-LDSC_reference_files

# conda
source $tool_dir/anaconda3/bin/activate
conda activate ldsc

# cluster
cells=(bigCTRL HSC GP Granulocyte MEMP-t MEMP MEP MEMP-Mast-Ery MEMP-Ery Early-Ery Late-Ery MEMP-MK MK MastP-t MastP Mast MDP Monocyte Kupffer cDC1 cDC2 pDC ASDC LMPP LP Cycling-LP PreProB ProB-1 ProB-2 Large-PreB Small-PreB IM-B NK ILCP T Hepatocyte Endothelia)

for d in FDR5_LFC1_padding_1k
do
  for cell in "${cells[@]}"
  do
    for chr in {1..22}
    do
      # Creating an annot file
      python $path2ldsc/make_annot.py \
        --bimfile $ldsc_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr}.bim \
        --bed-file $peak_dir/$d/$cell/${cell}_peak.hg19.bed \
        --annot-file $peak_dir/$d/$cell/${cell}.${chr}.annot.gz

      # Computing LD scores with an annot file
      python $path2ldsc/ldsc.py \
        --l2 \
        --bfile $ldsc_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr} \
        --ld-wind-cm 1 \
        --thin-annot \
        --annot $peak_dir/$d/$cell/$cell.${chr}.annot.gz \
        --print-snps $ldsc_files/hm3_no_MHC.list.txt \
        --out $peak_dir/$d/$cell/$cell.${chr}
    done
  done
done

