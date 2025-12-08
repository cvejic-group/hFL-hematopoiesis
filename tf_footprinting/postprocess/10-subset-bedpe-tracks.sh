#!/bin/bash

GOI="HLA-DRB1"
TSS="32589848"

BEDPE_IN="/work/DevM_analysis/data/E2G/glmmPG/HSC.glmmPG_FDR20.bedpe"
BEDPE_OUT="/work/aaa/projects/chrombpnet-devmult/scratch/bedpe/HSC.glmmPG_FDR20.${GOI}_tss.bedpe"


grep ${GOI} ${BEDPE_IN} |\
  awk -v tss="$TSS" 'BEGIN{OFS="\t"} {$5=tss; $6=tss; print}' \
  > ${BEDPE_OUT}
