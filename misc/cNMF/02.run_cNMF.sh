#!/usr/bin/env bash

# nohup bash 02.run_cNMF.sh >> logs/02.run_cNMF.log &

# set up
work_dir=$HOME/project/20231127_DevM/cNMF
data_dir=$work_dir/data
plot_dir=$work_dir/plots
cell=FL_HSC

# conda
source ~/software/anaconda3/bin/activate
conda activate cNMF

input=$data_dir/${cell}_batchcorrect

cnmf prepare --seed 14 --n-iter 100 --seed 14 -k `seq 2 30` --output-dir $data_dir --name $cell -c $input.Corrected.HVG.Varnorm.h5ad --genes-file $input.Corrected.HVGs.txt --tpm $input.TP10K.h5ad

# limit BLAS thread usage for cnmf thread
# otherwise each cnmf will use 8 threads (CPU cores)
OMP_NUM_THREADS=8 parallel -j 4 cnmf factorize --output-dir $data_dir --name $cell --worker-index {} ::: {0..3}

# --skip-completed-runs not working
OMP_NUM_THREADS=16 cnmf combine --output-dir $data_dir --name $cell

# k selection plot
OMP_NUM_THREADS=16 cnmf k_selection_plot --output-dir $data_dir --name $cell

# consensus
# k 5
OMP_NUM_THREADS=24 cnmf consensus --output-dir $data_dir --name $cell --components 5 --local-density-threshold 0.01 --show-clustering --build-reference

# k 7
OMP_NUM_THREADS=24 cnmf consensus --output-dir $data_dir --name $cell --components 7 --local-density-threshold 0.01 --show-clustering --build-reference

# k 12
OMP_NUM_THREADS=24 cnmf consensus --output-dir $data_dir --name $cell --components 12 --local-density-threshold 0.02 --show-clustering --build-reference
