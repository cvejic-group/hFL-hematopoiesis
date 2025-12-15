# Requirements:
# Environment should contain: https://github.com/aertslab/pycisTopic/tree/polars_1xx

binary_matrix_path='~/work/MalletRun/tmp/binary_accessibility_matrix.mtx'
outfile='~/work/MalletRun/tmp/corpus'
mallet_path='~/work/Mallet/bin/mallet'

pycistopic topic_modeling mallet create_corpus -i {binary_matrix_path} -o {outfile} -m 700 -b {mallet_path}
