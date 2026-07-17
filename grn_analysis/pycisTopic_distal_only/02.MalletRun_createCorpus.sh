#!/usr/bin/env bash
# Step 2: create the Mallet corpus from the binary accessibility matrix.
# Requirements: environment with pycisTopic polars_1xx
#   https://github.com/aertslab/pycisTopic/tree/polars_1xx
set -euo pipefail

# ---- paths (must match 01.build_cistopic_distal.py) ----
RUN_DIR='~/Mallet_res'
mallet_path='/work/Local_Data/utilities/Mallet/bin/mallet'

binary_matrix_path="${RUN_DIR}/tmp/binary_accessibility_matrix.mtx"
outfile="${RUN_DIR}/tmp/corpus"

# -m = memory in GB
pycistopic topic_modeling mallet create_corpus \
    -i "${binary_matrix_path}" \
    -o "${outfile}" \
    -m 700 \
    -b "${mallet_path}"
