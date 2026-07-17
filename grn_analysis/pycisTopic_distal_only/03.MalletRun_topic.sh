#!/usr/bin/env bash
# Step 3: run Mallet topic modeling over a range of topic numbers.
set -euo pipefail

RUN_DIR='/users/chenweiy/local_proj/Dev_M/00.Revision1_20260418/20260518_UMAP_with_distal_peaks/results/01.FL_distal_pycisTopic/Mallet_res'
RUN_NAME='01.FL_distal_pycisTopic'
mallet_path='/work/Local_Data/utilities/Mallet/bin/mallet'

mallet_corpus_path="${RUN_DIR}/tmp/corpus.mallet"
output_prefix="${RUN_DIR}/${RUN_NAME}"

# -t topic numbers | -p cpu | -n iterations | -a alpha | -A alpha_by_topic
# -e eta | -E eta_by_topic | -s random_state | -m memory(GB) | -v verbose
pycistopic topic_modeling mallet run \
    -i "${mallet_corpus_path}" \
    -o "${output_prefix}" \
    -t 15 20 25 30 35 40 45 \
    -p 30 \
    -n 150 \
    -a 50 -A True \
    -e 0.1 -E False \
    -s 555 \
    -m 1200 \
    -b "${mallet_path}" \
    -v
