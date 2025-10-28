#!/usr/bin/env Python

# u3-gpu-2
# conda activate scvi-py312
# python 03.multivi.py
# nohup python 03.multivi.py > logs/03.multivi.log &

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

work_dir = "/work/DevM_analysis/01.annotation/11.subclustering/BlineageByPCW"
data_dir = f"{work_dir}/data"
plot_dir = f"{work_dir}/plots"


for pcw in range(5, 19):
    print(f"PCW{pcw}")
    # load data
    print("Merge RNA/ATAC for multiVI")
    mdata = mu.read(f"{data_dir}/PCW{pcw}_wnn.h5mu")
    rna = mdata["rna"][:, mdata["rna"].var["highly_variable"]].copy()
    atac = mdata["atac"][:, mdata["atac"].var["selected"]].copy()
    atac.layers["counts"] = atac.X.copy()
    adata_paired = ad.concat([rna.copy(), atac.copy()], axis="var")
    adata_paired.obs = adata_paired.obs.join(atac.obs[["donorID", "libraryID"]])
    adata_paired.obs["modality"] = "paired"
    adata_paired

    # save
    adata_paired.write_h5ad(f"{data_dir}/PCW{pcw}_adata_paired.h5ad")

    print("Prepare objects")
    adata_mvi = scvi.data.organize_multiome_anndatas(adata_paired)

    print("Set up multiVI")
    if (pcw in ["5", "6", "7", "14", "17"]):
        scvi.model.MULTIVI.setup_anndata(
            adata_mvi,
            batch_key="modality",
            categorical_covariate_keys = ["libraryID", "donorID"],
            layer = "counts",
        )
    else:
        scvi.model.MULTIVI.setup_anndata(
            adata_mvi,
            batch_key="modality",
            categorical_covariate_keys = ["libraryID"],
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
    mvi.save(f"{data_dir}/PCW{pcw}_model_multiVI", overwrite=True, save_anndata=True)

    # get multiVI
    adata = mdata["rna"]
    SCVI_LATENT_KEY = "X_multiVI"
    adata.obsm[SCVI_LATENT_KEY] = mvi.get_latent_representation()
    sc.pp.neighbors(adata, use_rep = "X_multiVI")
    sc.tl.umap(adata)

    # save
    adata.write_h5ad(f"{data_dir}/PCW{pcw}_multiVI.h5ad")
    adata.obs.to_csv(f"{data_dir}/PCW{pcw}_mvi_cellmeta.csv", index=True)

    # umap
    umap_df = pd.DataFrame(adata.obsm['X_umap'], index=adata.obs_names, columns=['mvi_UMAP_1', 'mvi_UMAP_2'])
    umap_df = umap_df.rename_axis('barcode').reset_index()
    umap_df.to_csv(f"{data_dir}/PCW{pcw}_mvi_umap.csv")


