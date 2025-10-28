#!/usr/bin/env python

'''
conda activate scanpy
nohup python 00_FL_wnn.py > 00_FL_wnn.log &
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
        list_leiden_res = [0.5, 1, 2, 3, 4, 5, 6]
        print("Leiden clustering", file=sys.stderr)
        for r in list_leiden_res:
            print(r)
            # use igraph and fixed iterations to speed up
            sc.tl.leiden(mdata, resolution=r, key_added=f'leiden_wnn_{r}'.format(r))
    if run_diffmap:
        print("Diffusion map", file=sys.stderr)
        sc.tl.diffmap(mdata)
    return mdata


def run_dge(mdata=None, group=None, outdir=None):
    rna = mdata['rna']
    rna.obs[group] = mdata.obs[group].copy()
    sc.tl.rank_genes_groups(rna, groupby=group, method="wilcoxon")
    marker_df = sc.get.rank_genes_groups_df(rna, group=None)
    marker_df.to_csv(f"{outdir}/data/FL_wnn_markerGenes.v00.csv", index=False)
    # top 50
    result = rna.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    marker_df_top = pd.DataFrame({'cluster_' + group: result[key][group] for group in groups for key in ['names']}).head(50)
    marker_df_top.to_csv(f"{outdir}/data/FL_wnn_markerGenes_top50.v00.csv", index=False)
    return mdata


# variables
work_dir = '/work/home/project/20231127_DevM/multiome_wnn'
run = 'multiome_48FL'

# outdir
outdir = f"{work_dir}/{run}"
if not os.path.exists(outdir):
    os.makedirs(outdir)
if not os.path.exists(f"{outdir}/data/"):
    os.makedirs(f"{outdir}/data/")
if not os.path.exists(f"{outdir}/plots/"):
    os.makedirs(f"{outdir}/plots/")

# set up
rna = sc.read_h5ad(
    "/work/home/project/20231127_DevM/devm_rproj/output/intg_rna_48FL.00/rna_clustered.v00.h5ad"
)
print(rna, file=sys.stderr)

atac = sc.read_h5ad(
    "/work/home/project/20231127_DevM/snapatac/atac_clustered.v00.h5ad"
)
print(atac, file=sys.stderr)

# build
mdata = build_mdata(rna=rna, atac=atac)
mdata = run_wnn(mdata=mdata, L2norm=False)
mdata = post_dimred_processing(mdata=mdata)
mdata = run_dge(mdata=mdata, group='leiden_wnn_3', outdir=outdir)
mdata.write(f"{outdir}/data/FL_wnn_clustered.v00.h5mu")
mdata.obs.to_csv(f"{outdir}/data/FL_wnn_cellmeta.v00.csv", index=True)

print("Preprocessing completed!", file=sys.stderr)


