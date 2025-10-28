#!/usr/bin/env python

import gc
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.backends.backend_pdf as mpdf
from matplotlib.pyplot import rc_context

import anndata as ad
import scanpy as sc
import muon as mu

work_dir = '/work/DevM_analysis/02.abundance/Milo_FL_PCW250401'
data_dir = f"{work_dir}/data"
dataset = "FL_wnn"
mdata_raw = mu.read(
    f"{data_dir}/{dataset}_clustered.h5mu"
)
mdata_raw

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

k = 100
adata = prepare_data(mdata = mdata_raw, k = k)
adata.write('data/adata_4_milo.h5ad')


