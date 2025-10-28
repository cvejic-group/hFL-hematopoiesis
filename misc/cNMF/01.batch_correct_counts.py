#!/usr/bin/env python

'''
purpose: do batch correction on counts for cNMF

usage:
conda activate cNMF
nohup python 01.batch_correct_counts.py > logs/01.batch_correct_counts.log &
'''

import scanpy as sc
from cnmf import cNMF, Preprocess

# set up
work_dir = "/work/home/project/20231127_DevM/cNMF"
data_dir = f"{work_dir}/data"
plot_dir = f"{work_dir}/plots"
cell = "FL_HSC"

# load
adata = sc.read_h5ad("/work/DevM_analysis/01.annotation/11.subclustering/HSC/data/FL_rna_clustered.v01.h5ad")
print(adata)

#Initialize the Preprocess object
p = Preprocess(random_seed=14)

#Batch correct the data and save the corrected high-variance gene data to adata_c, and the TPM normalized data to adata_tpm
(adata_c, adata_tpm, hvgs) = p.preprocess_for_cnmf(adata, harmony_vars=['donorID', 'libraryID'], n_top_rna_genes = 2000,
                                                   librarysize_targetsum= 1e6,
                                                   save_output_base=f'{data_dir}/{cell}_batchcorrect')


