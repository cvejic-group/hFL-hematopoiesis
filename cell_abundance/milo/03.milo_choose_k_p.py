#!/usr/bin/env python

# conda activate pertpy
# nohup python 03.milo_choose_k_p.py > 03.milo_choose_k_p.log &

import warnings
from numba.core.errors import NumbaDeprecationWarning
warnings.filterwarnings(action='once')
warnings.simplefilter(action='once')
warnings.simplefilter(action="ignore", category=NumbaDeprecationWarning)
warnings.simplefilter(action="ignore", category=FutureWarning)
warnings.simplefilter(action="ignore", category=DeprecationWarning)

import os
import gc
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as mpdf
from matplotlib.pyplot import rc_context

import anndata as ad
import scanpy as sc
import muon as mu
import pertpy as pt
milo = pt.tl.Milo()

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=100, frameon=False, figsize=(8, 7), facecolor="white")
sc.logging.print_versions()

# dirs
work_dir = '/work/DevM_analysis/02.abundance/Milo_FL_PCW250401'
data_dir = f"{work_dir}/data"
dataset = "FL_wnn"
new_anno = "anno_wnn_v51"

# load
mdata = mu.read(
    f"{data_dir}/{dataset}_clustered.h5mu"
)
mdata

def prepare_data(mdata = None, k = None):
    mu.pp.neighbors(mdata, n_neighbors = k)
    # fake adata
    adata = ad.AnnData(mdata['rna'].X)
    # obs
    adata.obs = mdata.obs.copy()
    # graph
    adata.uns = mdata.uns.copy()
    adata.obsp = mdata.obsp.copy()
    # umap
    adata.obsm['X_umap'] = mdata.obsm['X_umap'].copy()
    # harmony from atac (4 graph refinement in milo)
    adata.obsm['X_harmony'] = mdata['atac'].obsm['X_harmony'].copy()
    adata.uns["neighbors"]["params"]["use_rep"] = 'X_harmony'
    return adata

def optimize_k_p(adata=None, k = None, p = None):
    # make_nhoods
    milo.make_nhoods(adata, prop=p)
    nhood_size = adata.obsm["nhoods"].toarray().sum(0)
    median, mean = str(round(np.median(nhood_size), 1)), str(round(np.mean(nhood_size), 1))
    mini, maxi = str(round(np.min(nhood_size), 1)), str(round(np.max(nhood_size), 1))
    #fig = plt.hist(nhood_size, bins=50)
    ax = sns.histplot(nhood_size, bins=50)
    ax.set(xlabel='# cells in neighbourhood', ylabel='# neighbouthoods', title = f"k={k}, prop={p}")
    ax.text(x = 0.95, y = 0.95, s = f"mean: {mean}", horizontalalignment='right', transform=ax.transAxes)
    ax.text(x = 0.95, y = 0.85, s = f"median: {median}", horizontalalignment='right', transform=ax.transAxes)
    ax.text(x = 0.95, y = 0.75, s = f"range: [{mini}, {maxi}]", horizontalalignment='right', transform=ax.transAxes)
    return ax

# p lst
p_lst = [0.05, 0.1, 0.2, 0.3]
# k lst
k_lst = [30, 50, 75, 100, 150, 200]

# for
for k in k_lst:
    adata = prepare_data(mdata, k = k)
    for p in p_lst:
        ax = optimize_k_p(adata, k = k, p = p)
        plt.savefig(f"plots/Milo_nh_dist.k={k}.p={p}.pdf", format="pdf", bbox_inches="tight")
        plt.close()

