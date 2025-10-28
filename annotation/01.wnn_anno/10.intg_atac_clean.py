#!/usr/bin/env python

# usage:
# on JupyterLab 418
# conda activate scanpy
# nohup python 01.intg_atac_clean.py > 01.intg_atac_clean.log &

import os
import sys
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
work_dir = "/work/DevM_analysis/01.annotation/10.integration_joint_clean"
data_dir = f"{work_dir}/data"
plot_dir = f"{work_dir}/plots"
atac = data.subset(obs_indices=pd.Series(data.obs_names)[data.obs[new_anno]!="LowQ"], out=f"{data_dir}")
atac = atac[0]
atac.obsm = None # slim

# close old
data.close()

# spectral
snap.pp.select_features(atac, n_features=50000)
print("spectral")
snap.tl.spectral(atac) #Default: n_comps=30, weighted_by_sd=True 

# batch
print("harmony")
adata = atac.to_adata()
atac.obsm["X_harmony"] = harmonize(adata.obsm["X_spectral"][:, :],
                                   adata.obs, batch_key=['libraryID', 'donorID'])
atac.obsm["X_harmony"] = np.float64(atac.obsm["X_harmony"])

# umap
print("umap")
snap.tl.umap(atac, use_rep="X_harmony")

# cluster
print("leiden")
snap.pp.knn(atac, use_rep="X_harmony")
list_leiden_res = [0.5, 1, 2, 3, 4, 5, 6]
for r in list_leiden_res:
    # iter = 2 to speed up
    snap.tl.leiden(atac, resolution=r, key_added=f'leiden_atac_{r}'.format(r), n_iterations=2)

# save
adata = atac.to_adata()

# level
adata.obs[new_anno] = adata.obs[new_anno].astype("category")
anno_order = ["HSC", "35-MDP?", "GP", "Granulocyte",
              "49-MEMP", "MEMP", "MEP", "MEMP-Mast-Ery", "MEMP-Ery", "Early-Ery", "Late-Ery",
              "MEMP-MK", "MK", "49:2-MastP", "MastP", "Mast",
              "Monocyte", "Kupffer", "cDC1", "cDC2", "pDC", "ASDC",
              "LMPP", "LP", "Cycling_LP", "PreProB", "ProB-1", "ProB-2", "Large-PreB", "Small-PreB", "IM-B",
              "NK", "ILCP", "T",
              "Hepatocyte", "Endothelia"]
adata.obs[new_anno] = adata.obs[new_anno].cat.reorder_categories(anno_order)

print("saving")
adata.write_h5ad(f"{data_dir}/FL_atac_snapatac2.h5ad")
atac.close()




