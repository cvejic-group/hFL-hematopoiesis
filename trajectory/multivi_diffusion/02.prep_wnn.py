#!/usr/bin/env python

'''
on JupyterLab 418
conda activate muon
nohup python 02.prep_wnn.py > logs/02.prep_wnn.log &
'''

import os
import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import scanpy as sc
import muon as mu
import warnings
warnings.filterwarnings(action='once')
warnings.simplefilter(action='once')

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=100, frameon=False, figsize=(8, 7), facecolor="white")
sc.logging.print_versions()


def build_mdata(rna=None, atac=None):
    # reorder rna (since atac is huge)
    rna = rna[atac.obs_names].copy()
    assert (rna.obs_names == atac.obs_names).all()
    # normalize RNA (for plotting purpose) --  h5ads were from seurat, adata.X is counts
    rna.layers["counts"] = rna.X.copy()
    sc.pp.normalize_total(rna, target_sum=1e4)
    sc.pp.log1p(rna)
    # mudata
    mdata = mu.MuData({"rna": rna, "atac": atac})
    return(mdata)


def run_wnn(mdata=None, use_rep='X_harmony', L2norm=False):
    if L2norm:
        mu.pp.l2norm(mdata['rna'], rep=use_rep)
        mu.pp.l2norm(mdata['atac'], rep=use_rep)
    # wnn
    sc.pp.neighbors(mdata['rna'], use_rep=use_rep)
    sc.pp.neighbors(mdata['atac'], use_rep=use_rep)
    mu.pp.neighbors(mdata)
    return(mdata)


def post_dimred_processing(mdata=None, run_clustering=True, run_diffmap=False):
    print("UMAP", file=sys.stderr)
    mu.tl.umap(mdata, random_state=10)
    # Leiden clustering
    if run_clustering:
        list_leiden_res = [0.1, 0.5, 1, 1.5, 2, 3, 4]
        print("Leiden clustering", file=sys.stderr)
        for r in list_leiden_res:
            print(r)
            # use igraph and fixed iterations to speed up
            sc.tl.leiden(mdata, resolution=r, key_added=f'leiden_wnn_{r}'.format(r), flavor="igraph", n_iterations=2)
    if run_diffmap:
        print("Diffusion map", file=sys.stderr)
        sc.tl.diffmap(mdata)
    return mdata


def run():
    # work dir
    work_dir = "/work/DevM_analysis/01.annotation/11.subclustering/BlineageByPCW"
    data_dir = f"{work_dir}/data"
    plot_dir = f"{work_dir}/plots"

    # pcw
    for pcw in [str(i) for i in range(5, 19)]:
        print(pcw)
        # set up
        rna = sc.read_h5ad(
            f"{data_dir}/PCW{pcw}_rna.h5ad"
        )
        print(rna, file=sys.stderr)

        atac = sc.read_h5ad(
            f"{data_dir}/PCW{pcw}_atac.h5ad"
        )
        print(atac, file=sys.stderr)

        # build
        mdata = build_mdata(rna=rna, atac=atac)
        mdata = run_wnn(mdata=mdata)
        mdata = post_dimred_processing(mdata=mdata, run_clustering=False)
        # obs
        mdata.obs['libraryID'] = mdata['atac'].obs['libraryID'].copy()
        mdata.obs['donorID'] = mdata['atac'].obs['donorID'].copy()
        mdata.obs['sampleID'] = mdata['atac'].obs['sampleID'].copy()
        mdata.obs['Sex'] = mdata['rna'].obs['Sex'].copy()
        mdata.obs['PCW'] = mdata['atac'].obs['PCW'].copy()
        mdata.obs['Batch'] = mdata['atac'].obs['Batch'].copy()
        mdata.obs['anno_wnn_v51'] = mdata['atac'].obs['anno_wnn_v51'].copy()
        # rna
        rna = mdata['rna']
        rna.obs['libraryID'] = mdata.obs['libraryID'].copy()
        rna.obs['donorID'] = mdata.obs['donorID'].copy()
        rna.obs['sampleID'] = mdata.obs['sampleID'].copy()
        rna.obs['PCW'] = mdata.obs['PCW'].copy()
        rna.obs['Batch'] = mdata.obs['Batch'].copy()
        rna.obs['anno_wnn_v51'] = mdata.obs['anno_wnn_v51'].copy()
        # save
        mdata.write(f"{data_dir}/PCW{pcw}_wnn.h5mu")
        mdata.obs.to_csv(f"{data_dir}/PCW{pcw}_wnn_cellmeta.csv", index=True)
        # umap
        umap_df = pd.DataFrame(mdata.obsm['X_umap'], index=mdata.obs_names, columns=['WNN_UMAP_1', 'WNN_UMAP_2'])
        umap_df = umap_df.rename_axis('barcode').reset_index()
        umap_df.to_csv(f"{data_dir}/PCW{pcw}_wnn_umap.csv")


if __name__ == "__main__":
    try:
        run()
        print("Complete!", file=sys.stderr)
    except KeyboardInterrupt:
        print("User interrupted me! ;-) Bye!")
        sys.exit(0)


