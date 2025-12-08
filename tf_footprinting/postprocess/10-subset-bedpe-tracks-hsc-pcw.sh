#!/bin/bash

GOI="CDCA3"

for pcw in {5..18}; do
    echo "HSC_PCW${pcw}-${GOI}"
    BEDPE_IN="/work/aaa/projects/chrombpnet-devmult/pipeline/results/seqlet-target/glmmPG/hsc_unified/HSC_PCW${pcw}_seqlet_glmmPG.counts_scores.bedpe.gz"
    BEDPE_OUT="/work/aaa/projects/chrombpnet-devmult/scratch/bedpe/HSC_PCW${pcw}_seqlet_glmmPG.counts_scores.${GOI}.bedpe"
    gzip -cd ${BEDPE_IN} | awk "\$18 ~ /${GOI}/" > ${BEDPE_OUT}
done
