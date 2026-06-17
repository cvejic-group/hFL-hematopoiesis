#!/usr/bin/env Python


# u3-gpu-4
# conda activate scvi-py312
# python 04.run_multivi.py

import sys, os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import anndata as ad
import scanpy as sc
import muon as mu
import scvi
import torch
torch.set_float32_matmul_precision("high")

scvi.settings.seed = 0
print("scvi-tools version:", scvi.__version__)

work_dir = "/work/DevM_analysis/01.annotation/11.subclustering"
celltype = "HSC"
data_dir = f"{work_dir}/{celltype}/data"
plot_dir = f"{work_dir}/{celltype}/plot"

print("Merge RNA/ATAC for multiVI")
mdata = mu.read(f"{data_dir}/FL_wnn_clustered.v00.h5mu")
rna = mdata["rna"][:, mdata["rna"].var["highly_variable"]].copy()
atac = mdata["atac"][:, mdata["atac"].var["selected"]].copy()
atac.layers["counts"] = atac.X.copy()
adata_paired = ad.concat([rna.copy(), atac.copy()], axis="var")
adata_paired.obs = adata_paired.obs.join(atac.obs[["donorID", "libraryID"]])
adata_paired.obs["modality"] = "paired"
adata_paired

# save
adata_paired.write_h5ad(f"{data_dir}/FL_adata_paired.h5ad")

print("Prepare objects")
adata_mvi = scvi.data.organize_multiome_anndatas(adata_paired)

print("Set up multiVI")
scvi.model.MULTIVI.setup_anndata(
    adata_mvi,
    batch_key="modality", 
    categorical_covariate_keys = ["libraryID", "donorID"],
    layer = "counts",
)
mvi = scvi.model.MULTIVI(
    adata_mvi,
    n_genes = len(rna.var_names),
    n_regions = len(atac.var_names),
)
mvi.view_anndata_setup()

print("Train multiVI")
mvi.train()

print("Save")
mvi.save(f"{data_dir}/model_multiVI", overwrite=True, save_anndata=True)


# get multiVI
adata = mdata["rna"]
SCVI_LATENT_KEY = "X_multiVI"
adata.obsm[SCVI_LATENT_KEY] = mvi.get_latent_representation()
sc.pp.neighbors(adata, use_rep = "X_multiVI")
sc.tl.umap(adata)

# leiden
list_leiden_res = [0.1, 0.3, 0.5, 1, 1.5, 2, 3, 4]
print("Leiden clustering", file=sys.stderr)
for r in list_leiden_res:
    print(r)
    # could use igraph and fixed iterations to speed up
    sc.tl.leiden(adata, resolution=r, key_added=f'leiden_mvi_{r}'.format(r))

# save
adata.write_h5ad(f"{data_dir}/FL_multiVI.h5ad")
adata.obs.to_csv(f"data/FL_multiVI_cellmeta.csv", index=True)

