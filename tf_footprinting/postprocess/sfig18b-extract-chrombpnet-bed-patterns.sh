DIR_BASE="/work/aaa/projects/chrombpnet-devmult/"
DIR_RESULTS="${DIR_BASE}/pipeline/results/chrombpnet_nobias/pretrained_bias"
DIR_OUT="${DIR_BASE}/scratch/bed"

## EMOMES
CELLTYPE="NK"
TF_NAME="EOMES"
PATTERN_TYPE="pos"
PATTERN_IDX="4"

PATTERN_REGEX="${PATTERN_TYPE}_patterns\.pattern_${PATTERN_IDX}\."
BED_IN="${DIR_RESULTS}/${CELLTYPE}/mean/modisco/${CELLTYPE}_mean.counts_scores.sorted.bed.gz"
BED_OUT="${DIR_OUT}/chrombpnet-${CELLTYPE}-${PATTERN_TYPE}-${PATTERN_IDX}-${TF_NAME}.bed"

gzip -cd ${BED_IN} | awk "\$4 ~ /${PATTERN_REGEX}/" > ${BED_OUT}


# CEBPx (D,E,G,B,A)
CELLTYPE="Monocyte"
TF_NAME="CEBPx"
PATTERN_TYPE="pos"
PATTERN_IDX="1"

PATTERN_REGEX="${PATTERN_TYPE}_patterns\.pattern_${PATTERN_IDX}\."
BED_IN="${DIR_RESULTS}/${CELLTYPE}/mean/modisco/${CELLTYPE}_mean.counts_scores.sorted.bed.gz"
BED_OUT="${DIR_OUT}/chrombpnet-${CELLTYPE}-${PATTERN_TYPE}-${PATTERN_IDX}-${TF_NAME}.bed"

gzip -cd ${BED_IN} | awk "\$4 ~ /${PATTERN_REGEX}/" > ${BED_OUT}


# EBF1
CELLTYPE="IM-B"
TF_NAME="EBF1"
PATTERN_TYPE="pos"
PATTERN_IDX="12"

PATTERN_REGEX="${PATTERN_TYPE}_patterns\.pattern_${PATTERN_IDX}\."
BED_IN="${DIR_RESULTS}/${CELLTYPE}/mean/modisco/${CELLTYPE}_mean.counts_scores.sorted.bed.gz"
BED_OUT="${DIR_OUT}/chrombpnet-${CELLTYPE}-${PATTERN_TYPE}-${PATTERN_IDX}-${TF_NAME}.bed"

gzip -cd ${BED_IN} | awk "\$4 ~ /${PATTERN_REGEX}/" > ${BED_OUT}


# GATA
CELLTYPE="Late-Ery"
TF_NAME="GATA"
PATTERN_TYPE="pos"
PATTERN_IDX="0"

PATTERN_REGEX="${PATTERN_TYPE}_patterns\.pattern_${PATTERN_IDX}\."
BED_IN="${DIR_RESULTS}/${CELLTYPE}/mean/modisco/${CELLTYPE}_mean.counts_scores.sorted.bed.gz"
BED_OUT="${DIR_OUT}/chrombpnet-${CELLTYPE}-${PATTERN_TYPE}-${PATTERN_IDX}-${TF_NAME}.bed"

gzip -cd ${BED_IN} | awk "\$4 ~ /${PATTERN_REGEX}/" > ${BED_OUT}
