#!/usr/bin/env bash

# nohup bash 02.peak_liftOver.sh > 02.peak_liftOver.log &

work_dir=~/project/20231127_DevM/gwas/sLDSC
peak_dir=$work_dir/atac_peaks

d=FDR5_LFC1_padding_1k

for cell in CTRL HSC GP Granulocyte MEMP-t MEMP MEP MEMP-Mast-Ery MEMP-Ery Early-Ery Late-Ery MEMP-MK MK MastP-t MastP Mast MDP Monocyte Kupffer cDC1 cDC2 pDC ASDC LMPP LP Cycling-LP PreProB ProB-1 ProB-2 Large-PreB Small-PreB IM-B NK ILCP T Hepatocyte Endothelia
do
    # liftover
    liftOver $peak_dir/$d/$cell/${cell}_peak.hg38.bed ~/RefData/ucsc/liftOver/hg38ToHg19.over.chain.gz $peak_dir/$d/$cell/${cell}_peak.hg38ToHg19.bed $peak_dir/$d/$cell/${cell}_peak.hg38ToHg19.unlifted
    # rm scaffolds
    # no need to remove chr prefix (sort for safe)
    # in make_annot.py, it will add chr prefix to the bim files (stupid design)
    awk '/^chr[0-9]*\t/' $peak_dir/$d/$cell/${cell}_peak.hg38ToHg19.bed | \
      sort -k1,1V -k2,2n > $peak_dir/$d/$cell/${cell}_peak.hg19.bed
done

# create a big CTRL set
# by combining peaks from Zhang_2021_adult_atlas and Domcke_2020_fetal_atlas
# as described in "Chromatin accessibility during human first-trimester neurodevelopment" paper
for d in FDR1_LFC1_padding_1k FDR5_LFC1_padding_1k
do
  cell=bigCTRL
  mkdir $peak_dir/$d/$cell
  cat $peak_dir/$d/CTRL/CTRL_peak.hg19.bed /work/DevM_analysis/data/refATAC/Domcke_2020_fetal_atlas/master_peak.bed3 /work/DevM_analysis/data/refATAC/Zhang_2021_adult_atlas/cCRE_hg19.bed | sort -k1,1V -k2,2n > $peak_dir/$d/$cell/${cell}_peak.hg19.bed
done

