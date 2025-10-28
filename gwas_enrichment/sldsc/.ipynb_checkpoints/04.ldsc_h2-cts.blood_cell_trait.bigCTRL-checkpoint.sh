#!/usr/bin/env bash

# constrain threads usage by ldsc
export OMP_NUM_THREADS=10
export OPENBLAS_NUM_THREADS=10
export MKL_NUM_THREADS=10
export NUMEXPR_NUM_THREADS=10

work_dir=~/project/20231127_DevM/gwas/sLDSC
peak_dir=$work_dir/atac_peaks
ss_dir=$work_dir/munge_sumstats
out_dir=$work_dir/ldsc_h2-cts_bigCTRL

# ldsc
tool_dir=$HOME/software
path2ldsc=$tool_dir/ldsc
ldsc_files=/work/DevM_analysis/data/S-LDSC_reference_files
BASELINE=$ldsc_files/baselineLD_v2.3/baselineLD.
WEIGHTS=$ldsc_files/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
FRQ=$ldsc_files/1000G_Phase3_frq/1000G.EUR.QC.

# conda
source $tool_dir/anaconda3/bin/activate
conda activate ldsc

# blood cell count
phenos=(pbmc_cd34 baso_ukb basop_ukb eo_ukb eop_ukb hct_ukb hgb_ukb hlr_ukb hlsrc_ukb irf_ukb lymph_ukb lymphp_ukb mch_ukb mchc_ukb mcv_ukb mono_ukb monop_ukb mpv_ukb mrv_ukb mscv_ukb neut_ukb neutp_ukb pct_ukb pdw_ukb plt_ukb rbc_ukb rdw_cv_ukb ret_ukb retf_ukb wbc_ukb)

for d in FDR1_LFC1_padding_1k FDR5_LFC1_padding_1k
do
  for pheno in "${phenos[@]}"
  do
    #pheno=baso_ukb
    echo $pheno
    sumstats=$ss_dir/$pheno.sumstats.gz
    python $path2ldsc/ldsc.py \
      --h2-cts $sumstats \
      --ref-ld-chr $BASELINE \
      --ref-ld-chr-cts $work_dir/$d.bigCTRL.ldcts \
      --w-ld-chr $WEIGHTS \
      --out $out_dir/$d/$pheno
  done

done
