#!/usr/bin/env python

'''
on JupyterLab 418
conda activate scanpy
nohup python 02.wnn_clustering.py > 02.wnn_clustering.log &
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


def run_wnn(mdata=None, rna_rep='X_harmony', atac_rep='X_harmony', L2norm=False):
    if L2norm:
        mu.pp.l2norm(mdata['rna'], rep=rna_rep)
        mu.pp.l2norm(mdata['atac'], rep=atac_rep)
    # wnn
    sc.pp.neighbors(mdata['rna'], use_rep=rna_rep)
    sc.pp.neighbors(mdata['atac'], use_rep=atac_rep)
    mu.pp.neighbors(mdata)
    return(mdata)


def post_dimred_processing(mdata=None, run_clustering=True, run_diffmap=False):
    print("UMAP", file=sys.stderr)
    mu.tl.umap(mdata, random_state=10)
    # Leiden clustering
    if run_clustering:
        list_leiden_res = [0.1, 0.3, 0.5, 1, 1.5, 2, 3, 4]
        print("Leiden clustering", file=sys.stderr)
        for r in list_leiden_res:
            print(r)
            # could use igraph and fixed iterations to speed up
            sc.tl.leiden(mdata, resolution=r, key_added=f'leiden_wnn_{r}'.format(r))
    if run_diffmap:
        print("Diffusion map", file=sys.stderr)
        sc.tl.diffmap(mdata)
    return mdata


def run_dge(mdata=None, group=None, outdir=None):
    rna = mdata['rna']
    sc.tl.rank_genes_groups(rna, groupby=group, method="wilcoxon")
    marker_df = sc.get.rank_genes_groups_df(rna, group=None)
    marker_df.to_csv(f"{outdir}/FL_wnn_markerGenes.v00.csv", index=False)
    # top 100
    result = rna.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    marker_df_top = pd.DataFrame({group: result[key][group] for group in groups for key in ['names']}).head(100)
    marker_df_top.to_csv(f"{outdir}/FL_wnn_markerGenes_top100.v00.csv", index=False)
    # top 500
    result = rna.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    marker_df_top = pd.DataFrame({group: result[key][group] for group in groups for key in ['names']}).head(500)
    marker_df_top.to_csv(f"{outdir}/FL_wnn_markerGenes_top100.v00.csv", index=False)
    return mdata


# variables
work_dir = "/work/DevM_analysis/01.annotation/10.integration_joint_clean"
data_dir = f"{work_dir}/data"
plot_dir = f"{work_dir}/plots"

# set up
rna = sc.read_h5ad(
    f"{data_dir}/FL_rna_clustered.v00.h5ad"
)
print(rna, file=sys.stderr)

atac = sc.read_h5ad(
    f"{data_dir}/FL_atac_snapatac2.h5ad"
)
print(atac, file=sys.stderr)

# build
mdata = build_mdata(rna=rna, atac=atac)
mdata = run_wnn(mdata=mdata)
mdata = post_dimred_processing(mdata=mdata, run_clustering=False)
mdata = run_dge(mdata=mdata, group='anno_wnn_v5', outdir=data_dir)

# some obs
mdata.obs['libraryID'] = mdata['rna'].obs['libraryID'].copy()
mdata.obs['donorID'] = mdata['rna'].obs['donorID'].copy()
mdata.obs['sampleID'] = mdata['rna'].obs['sampleID'].copy()
mdata.obs['PCW'] = mdata['rna'].obs['PCW'].astype(int).copy()
mdata.obs['anno_wnn_v5'] = mdata['rna'].obs['anno_wnn_v5'].copy()

# save
mdata.obs.to_csv(f"{data_dir}/FL_wnn_cellmeta.v00.csv", index=True)
mdata.write(f"{data_dir}/FL_wnn_clustered.v00.h5mu")

# umap wnn
umap_df = pd.DataFrame(mdata.obsm['X_umap'], index=mdata.obs_names, columns=['WNN_UMAP_1', 'WNN_UMAP_2'])
umap_df = umap_df.rename_axis('barcode').reset_index()
umap_df.to_csv(f"{data_dir}/FL_wnn_umap.csv", index=False)

# umap rna
umap_df = pd.DataFrame(mdata['rna'].obsm['X_umap'], index=mdata.obs_names, columns=['RNA_UMAP_1', 'RNA_UMAP_2'])
umap_df = umap_df.rename_axis('barcode').reset_index()
umap_df.to_csv(f"{data_dir}/FL_rna_umap.csv", index=False)

# umap atac
umap_df = pd.DataFrame(mdata['atac'].obsm['X_umap'], index=mdata.obs_names, columns=['ATAC_UMAP_1', 'ATAC_UMAP_2'])
umap_df = umap_df.rename_axis('barcode').reset_index()
umap_df.to_csv(f"{data_dir}/FL_atac_umap.csv", index=False)


print("Preprocessing completed!", file=sys.stderr)

