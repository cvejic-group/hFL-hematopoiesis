#!/usr/bin/env python

# usage:
# on JupyterLab 418
# conda activate scanpy
# nohup python 02.prep_atac.py > 02.prep_atac.log &

import os
import snapatac2 as snap
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from harmony import harmonize


# anno
new_anno = "anno_wnn_v5"
df_anno = pd.read_csv("/work/DevM_analysis/01.annotation/09.annotation_joint/data/FL_wnn_cellmeta.v05.csv")
df_anno.index = df_anno["Unnamed: 0"]

# dataset
data = snap.read_dataset("/work/DevM_analysis/01.annotation/07.integration_atac/data/FL_atac_snapatac2-integration.annDataset.h5ads")
data.obs[new_anno] = df_anno.loc[data.obs_names, new_anno]

# subset
work_dir = "/work/DevM_analysis/01.annotation/11.subclustering"
celltype = "HSC"
data_dir = f"{work_dir}/{celltype}/data"
plot_dir = f"{work_dir}/{celltype}/plots"
atac = data.subset(obs_indices=pd.Series(data.obs_names)[data.obs[new_anno]=="HSC"], out=f"{data_dir}")
atac = atac[0]
atac.obsm = None # slim

# spectral
snap.pp.select_features(atac, n_features=50000)
snap.tl.spectral(atac) #Default: n_comps=30, weighted_by_sd=True 

# batch
adata = atac.to_adata()
if (all(adata.obs['libraryID'] == adata.obs['donorID'])):
    atac.obsm["X_harmony"] = harmonize(adata.obsm["X_spectral"][:, :],
                                       adata.obs, batch_key=['libraryID'])
else:
    atac.obsm["X_harmony"] = harmonize(adata.obsm["X_spectral"][:, :],
                                       adata.obs, batch_key=['libraryID', 'donorID'])
atac.obsm["X_harmony"] = np.float64(atac.obsm["X_harmony"])

# umap
snap.tl.umap(atac, use_rep="X_harmony")

# cluster
snap.pp.knn(atac, use_rep="X_harmony")
list_leiden_res = [0.5, 1, 2, 3, 4, 5, 6]
for r in list_leiden_res:
    # iter = 2 to speed up
    snap.tl.leiden(atac, resolution=r, key_added=f'leiden_atac_{r}'.format(r), n_iterations=2)

# save
atac.obsm['fragment_paired'] = atac.adatas.obsm['fragment_paired'] # slow
adata = atac.to_adata()
adata.write_h5ad(f"{data_dir}/FL_atac_snapatac2.h5ad")
atac.close()




