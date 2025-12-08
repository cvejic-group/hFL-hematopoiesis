#!/bin/bash

mkdir -p frags
CSV_FRAGCT="frags/hsc_pcw_fragcts.csv"

echo "sample_id,n_frags" > ${CSV_FRAGCT}

for pcw in {5..18}; do
  SAMPLE_ID="HSC_PCW${pcw}"
  echo "${SAMPLE_ID}"
  BED="/work/aaa/projects/chrombpnet-devmult/pipeline/resources/peaks/${SAMPLE_ID}/${SAMPLE_ID}.no_blacklist.bed"
  FRAGS="/work/DevM_analysis/data/sinto/frag_per_cell/HSC/${SAMPLE_ID}.fragments.tsv.gz"
  CTS=$(bedtools intersect -a ${BED} -b ${FRAGS} -c | awk '{k += $11} END {print k}')
  echo "${SAMPLE_ID},${CTS}" >> ${CSV_FRAGCT}
done
