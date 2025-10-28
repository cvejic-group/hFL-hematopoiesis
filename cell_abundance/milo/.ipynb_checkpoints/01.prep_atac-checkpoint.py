#!/usr/bin/env python

# usage:
# conda activate muon
# nohup python 01.prep_atac.py > 01.prep_atac.log &

import os, sys
import snapatac2 as snap
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from harmony import harmonize


# input
adata = sc.read_h5ad('/work/DevM_analysis/01.annotation/11.subclustering/blood/data/FL_atac_snapatac2.h5ad')
print(adata, file=sys.stderr)

# dirs
work_dir = "/work/DevM_analysis/02.abundance/Milo_FL_PCW250401"
data_dir = f"{work_dir}/data"
plot_dir = f"{work_dir}/plots"

# subset
atac = adata[(~adata.obs['PCW'].isin(['Mixed']))].copy()
atac.obs['anno_wnn_v51'] = atac.obs['anno_wnn_v51'].cat.remove_unused_categories()
atac.obs['PCW'] = atac.obs['PCW'].cat.remove_unused_categories()
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
atac.obsm["X_harmony"] = harmonize(atac.obsm["X_spectral"][:, :],
                                   atac.obs, batch_key=['libraryID', 'donorID'])
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
atac.write_h5ad(f"{data_dir}/FL_atac_clustered.h5ad")

