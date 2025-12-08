source /work/aaa/.bashrc
conda activate deeptools

DIR_IN="/work/aaa/projects/chrombpnet-devmult/pipeline/results/chrombpnet_nobias/hsc_unified"
DIR_OUT="/work/aaa/projects/chrombpnet-devmult/scratch/bigwig"
CHROMS="/work/aaa/projects/chrombpnet-devmult/pipeline/resources/genome/hg38.chrom.sizes"

## EARLY: PCW 5 (copy only, no merging)
ID_GROUP="early"
BWS_IN="${DIR_IN}/HSC_PCW5/mean/preds/HSC_PCW5_mean._bias.bw"
BW_OUT="${DIR_OUT}/HSC_${ID_GROUP}_mean._bias.bw"

echo "Copying early..."
cp $BWS_IN $BW_OUT


## MIDDLE: PCW 6 - 10
ID_GROUP="middle"
BWS_IN=$(for i in {6..10}; do echo -n "${DIR_IN}/HSC_PCW${i}/mean/preds/HSC_PCW${i}_mean._bias.bw "; done)
BW_OUT="${DIR_OUT}/HSC_${ID_GROUP}_mean._bias.bw"

BG_SUM="${DIR_OUT}/tmp_pred_bias${ID_GROUP}_sum.bedGraph"
BG_MEAN="${DIR_OUT}/tmp_pred_bias${ID_GROUP}_mean.bedGraph"

echo "Summing ${ID_GROUP}..."
bigWigMerge -threshold=-1000.0 ${BWS_IN} ${BG_SUM}

echo "Averaging ${ID_GROUP}..."
awk '{{print $1, $2, $3, $4/5}}' ${BG_SUM} > ${BG_MEAN}
rm ${BG_SUM}

echo "Converting ${ID_GROUP} to bigWig..."
bedGraphToBigWig ${BG_MEAN} ${CHROMS} ${BW_OUT}
rm ${BG_MEAN}


# LATE: PCW 11 - 18
ID_GROUP="late"
BWS_IN=$(for i in {11..18}; do echo -n "${DIR_IN}/HSC_PCW${i}/mean/preds/HSC_PCW${i}_mean._bias.bw "; done)
BW_OUT="${DIR_OUT}/HSC_${ID_GROUP}_mean._bias.bw"

BG_SUM="${DIR_OUT}/tmp_pred_bias${ID_GROUP}_sum.bedGraph"
BG_MEAN="${DIR_OUT}/tmp_pred_bias${ID_GROUP}_mean.bedGraph"

echo "Summing ${ID_GROUP}..."
bigWigMerge -threshold=-1000.0 ${BWS_IN} ${BG_SUM}

echo "Averaging ${ID_GROUP}..."
awk '{{print $1, $2, $3, $4/8}}' ${BG_SUM} > ${BG_MEAN}
rm ${BG_SUM}

echo "Converting ${ID_GROUP} to bigWig..."
bedGraphToBigWig ${BG_MEAN} ${CHROMS} ${BW_OUT}
rm ${BG_MEAN}

echo "Done."
