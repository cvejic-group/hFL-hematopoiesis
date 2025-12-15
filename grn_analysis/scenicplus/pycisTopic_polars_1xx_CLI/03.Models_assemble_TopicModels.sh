binary_matrix_path='~/work/MalletRun/tmp/binary_accessibility_matrix.mtx'
cells_path='~/work/MalletRun/tmp/cellNames.txt'
regions_path='~/work/MalletRun/tmp/peakNames.txt'
output_prefix='~/work/MalletRun/MalletRun'

## assemble topic models
pycistopic topic_modeling mallet stats -i ${binary_matrix_path} -c ${cells_path} -r ${regions_path} -o ${output_prefix} -t 15 20 25 30 35 40 45 -v
