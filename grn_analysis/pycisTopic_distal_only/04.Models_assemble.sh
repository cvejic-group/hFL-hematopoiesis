#!/usr/bin/env bash
# Step 4: assemble the Mallet outputs into pycisTopic LDA model objects (stats).
# Produces one <prefix>.<n>_topics.model.pkl per topic number.
set -euo pipefail

RUN_DIR='~/Mallet_res'
RUN_NAME='01.FL_distal_pycisTopic'

binary_matrix_path="${RUN_DIR}/tmp/binary_accessibility_matrix.mtx"
cells_path="${RUN_DIR}/tmp/cellNames.txt"
regions_path="${RUN_DIR}/tmp/peakNames.txt"
output_prefix="${RUN_DIR}/${RUN_NAME}"

pycistopic topic_modeling mallet stats \
    -i "${binary_matrix_path}" \
    -c "${cells_path}" \
    -r "${regions_path}" \
    -o "${output_prefix}" \
    -t 15 20 25 30 35 40 45 \
    -v
