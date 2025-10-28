#!/usr/bin/env python

# usage:
# on JupyterLab 418
# conda activate snapatac2
# nohup python 01.prep_atac.py > logs/01.prep_atac.log &

import os, sys
import snapatac2 as snap
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from harmony import harmonize

# input
adata = sc.read_h5ad('/work/DevM_analysis/01.annotation/11.subclustering/Blineage/data/FL_atac_snapatac2.h5ad')
print(adata, file=sys.stderr)

# dirs
work_dir = "/work/DevM_analysis/01.annotation/11.subclustering"
new_anno = "anno_wnn_v51"
celltype = "BlineageByPCW"
data_dir = f"{work_dir}/{celltype}/data"
plot_dir = f"{work_dir}/{celltype}/plots"

# filter low
filter_low = 30

# subset
for pcw in [str(i) for i in range(5, 19)]:
    print(pcw)
    # subset pcw
    groups_to_include = [pcw]
    atac = adata[adata.obs['PCW'].isin(groups_to_include)].copy()
    # filter low
    cluster_counts = atac.obs[new_anno].value_counts()
    cluster_2_keep = cluster_counts[cluster_counts >= filter_low].index.tolist()
    atac = atac[atac.obs[new_anno].isin(cluster_2_keep)].copy()
    # remove unused category
    atac.obs['PCW'] = atac.obs['PCW'].cat.remove_unused_categories()
    atac.obs['anno_wnn_v51'] = atac.obs['anno_wnn_v51'].cat.remove_unused_categories()
    atac.obs['libraryID'] = atac.obs['libraryID'].cat.remove_unused_categories()
    atac.obs['donorID'] = atac.obs['donorID'].cat.remove_unused_categories()
    atac.obs['sampleID'] = atac.obs['sampleID'].cat.remove_unused_categories()
    atac.obs['Batch'] = atac.obs['Batch'].cat.remove_unused_categories()
    
    # slim
    atac.obsm = None
    print(atac, file=sys.stderr)
    
    # spectral
    snap.pp.select_features(atac, n_features=50000)
    snap.tl.spectral(atac) #Default: n_comps=30, weighted_by_sd=True 
    
    # batch
    if (pcw in ["5", "6", "7", "14", "17"]):
        atac.obsm["X_harmony"] = harmonize(atac.obsm["X_spectral"][:, :],
                                           atac.obs, batch_key=['donorID', 'libraryID'])
    else:
        atac.obsm["X_harmony"] = harmonize(atac.obsm["X_spectral"][:, :],
                                           atac.obs, batch_key=['libraryID'])
    atac.obsm["X_harmony"] = np.float64(atac.obsm["X_harmony"])
    
    # umap
    snap.tl.umap(atac, use_rep="X_harmony")
    
    # save
    atac.write_h5ad(f"{data_dir}/PCW{pcw}_atac.h5ad")

