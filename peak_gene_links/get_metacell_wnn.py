"""
single cell analysis helper functions for mudata with RNA and ATAC modalities
"""

import gc
import logging as log
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
import numpy as np
import pandas as pd
from sknetwork.hierarchy import Paris, LouvainHierarchy, cut_balanced
from scipy.spatial.distance import cosine

import scanpy as sc
import snapatac2 as snap
import muon as mu
from harmony import harmonize


# batch correction helper
def are_columns_not_equivalent(df: pd.DataFrame, col1: str, col2: str) -> bool:
    """
    Check whether two columns in a DataFrame are functionally equivalent,
    meaning each unique value in col1 maps to exactly one value in col2 and vice versa.
    """
    col1_to_col2 = df.groupby(col1)[col2].nunique().eq(1).all()
    col2_to_col1 = df.groupby(col2)[col1].nunique().eq(1).all()
    return not (col1_to_col2 and col2_to_col1)


def prep_rna(mdata=None):
    log.info(f"Prep RNA")
    rna = mdata['rna']
    # data were from WNN
    # normalize_total and log1p were done
    sc.pp.highly_variable_genes(rna)
    sc.pp.pca(rna, n_comps=30)
    # Add batch correction using Harmony
    covariates = []
    if len(rna.obs['libraryID'].unique()) > 1:
        covariates.append('libraryID')
        if len(rna.obs['donorID'].unique()) > 1 and are_columns_not_equivalent(rna.obs, 'libraryID', 'donorID'):
            covariates.append('donorID')
    if len(covariates) > 0:
        log.info(f"Harmony with covariates: {covariates}")
        rna.obsm["X_harmony"] = harmonize(rna.obsm["X_pca"][:, :],
                                          rna.obs, batch_key=covariates)
    else:
        log.info(f"No batch correction needed")
        rna.obsm["X_harmony"] = rna.obsm["X_pca"][:, :]
    rna.obsm["X_harmony"] = np.float64(rna.obsm["X_harmony"])
    return mdata


def prep_atac(mdata=None):
    log.info(f"Prep ATAC")
    atac = mdata['atac']
    # default: n_features=50000
    snap.pp.select_features(atac, n_features=50000)
    # default: n_comps=30, weighted_by_sd=True
    snap.tl.spectral(atac)
    # Add batch correction using Harmony
    covariates = []
    if len(atac.obs['libraryID'].unique()) > 1:
        covariates.append('libraryID')
        if len(atac.obs['donorID'].unique()) > 1 and are_columns_not_equivalent(atac.obs, 'libraryID', 'donorID'):
            covariates.append('donorID')
    if len(covariates) > 0:
        log.info(f"Harmony with covariates: {covariates}")
        atac.obsm["X_harmony"] = harmonize(atac.obsm["X_spectral"][:, :],
                                           atac.obs, batch_key=covariates)
    else:
        log.info(f"No batch correction needed")
        atac.obsm["X_harmony"] = atac.obsm["X_spectral"][:, :]
    atac.obsm["X_harmony"] = np.float64(atac.obsm["X_harmony"])
    return mdata


def prep_wnn(mdata=None, rna_rep='X_harmony', atac_rep='X_harmony', L2norm=False, **kwargs):
    log.info(f"Prep WNN")
    if L2norm:
        mu.pp.l2norm(mdata['rna'], rep=rna_rep)
        mu.pp.l2norm(mdata['atac'], rep=atac_rep)
    # wnn
    sc.pp.neighbors(mdata['rna'], use_rep=rna_rep, **kwargs)
    sc.pp.neighbors(mdata['atac'], use_rep=atac_rep, **kwargs)
    # if n_obs is smaller than default n_multineighbors = 200, error (https://github.com/scverse/muon/issues/86#issuecomment-2565219453)
    n_obs = mdata.n_obs
    if n_obs > 500:
        mu.pp.neighbors(mdata, **kwargs)
    elif n_obs < 500 & n_obs > 200:
        mu.pp.neighbors(mdata, n_neighbors=10, n_multineighbors=100, **kwargs)
    elif n_obs < 200 & n_obs > 100:
        mu.pp.neighbors(mdata, n_neighbors=10, n_multineighbors=50, **kwargs)
    elif n_obs < 100 & n_obs > 50:
        mu.pp.neighbors(mdata, n_neighbors=5, n_multineighbors=30, **kwargs)
    else:
        mu.pp.neighbors(mdata, n_neighbors=5, n_multineighbors=10, **kwargs)
    return(mdata)


def merge_small_clusters(labels=None, adjacency=None, min_cluster_size=None):
    """Merge small clusters based on graph connectivity (same principle as Paris()).
    """
    unique_labels, counts = np.unique(labels, return_counts=True)

    # Identify small clusters
    small_clusters = unique_labels[counts < min_cluster_size]

    if len(small_clusters) == 0:
        return labels  # No small clusters, return original labels

    # Compute cluster centroids based on graph adjacency
    cluster_centroids = {}
    for cluster in unique_labels:
        cluster_indices = np.where(labels == cluster)[0]
        cluster_centroids[cluster] = adjacency[cluster_indices].mean(axis=0)  # Mean adjacency row
    # Merge small clusters
    new_labels = labels.copy()
    for small_cluster in small_clusters:
        small_cluster_indices = np.where(labels == small_cluster)[0]
        # Compute connectivity (sum of adjacency values) to other clusters
        connectivity_scores = {
            cluster: np.sum(adjacency[small_cluster_indices][:, labels == cluster])
            for cluster in unique_labels if cluster not in small_clusters
        }
        # Find the most connected larger cluster
        if connectivity_scores:
            closest_cluster = max(connectivity_scores, key=connectivity_scores.get)
            new_labels[small_cluster_indices] = closest_cluster
    return new_labels


def get_metacell_groups(
    neighbor_matrix=None,
    max_group_size=10,
    min_group_size=5,
    method="Paris",
):
    """
    Divide observations by knn graph into clusters of a
    maximum size.
    Options for method: ['Paris', 'Louvain']
    """
    # check input size
    n_obs = neighbor_matrix.shape[0]
    if n_obs < max_group_size:
        warnings.warn(
            "n_obs is smaller than max_group_size: "
            f"{n_obs} < {max_group_size}. "
            f"Setting max_group_size to {n_obs}."
        )
        max_group_size = n_obs  # TODO: better choice?
    # create a dendrogram
    if method == "Paris":
        paris = Paris()
        dendrogram = paris.fit_predict(neighbor_matrix)
    elif method == "Louvain":
        louvain = LouvainHierarchy()
        dendrogram = louvain.fit_predict(neighbor_matrix)
    else:
        raise ValueError(f"Function get_metacell_groups: method '{method}' unknown")
    # cut into groups
    groups = cut_balanced(dendrogram, max_cluster_size = max_group_size)
    # add smaller clusters to closest cluster
    groups_final = merge_small_clusters(labels=groups, adjacency=neighbor_matrix, min_cluster_size=min_group_size)
    return groups_final


def get_metacells(mdata=None, run_wnn=False, **kwargs):
    # Error handling if "connectivities" doesn't exist
    if "connectivities" not in mdata.obsp:
        if run_wnn:
            mdata = prep_rna(mdata)
            mdata = prep_atac(mdata)
            mdata = prep_wnn(mdata)
        else:
            raise ValueError("No 'connectivities' matrix found in mdata.obsp")
    if mdata.obsp["connectivities"].shape[0] != mdata.n_obs:
        if run_wnn:
            mdata = prep_rna(mdata)
            mdata = prep_atac(mdata)
            mdata = prep_wnn(mdata)
        else:
            raise ValueError(f"Dimensions of connectivities matrix ({mdata.obsp['connectivities'].shape}) do not match number of observations ({mdata.n_obs})")
    n_mat = mdata.obsp["connectivities"].copy() # could be large
    groups = get_metacell_groups(n_mat, **kwargs)

    # Create dataframe with barcode and group columns
    df = pd.DataFrame({
        'barcode': mdata.obs.index,
        'mc_group': groups
    })
    return df


def get_metacells_by_group(mdata=None, group_by=None, groups=None, **kwargs):
    if group_by not in mdata.obs.columns:
        raise ValueError(f"group_by '{group_by}' not found in mdata.obs columns")

    # Check if libraryID exists in both modalities (for harmony)
    if 'libraryID' not in mdata['rna'].obs.columns:
        raise ValueError("'libraryID' not found in RNA modality")
    if 'libraryID' not in mdata['atac'].obs.columns:
        raise ValueError("'libraryID' not found in ATAC modality")

    # Check if donorID exists in both modalities (for harmony)
    if 'donorID' not in mdata['rna'].obs.columns:
        raise ValueError("'donorID' not found in RNA modality")
    if 'donorID' not in mdata['atac'].obs.columns:
        raise ValueError("'donorID' not found in ATAC modality")

    dfs = []
    if groups is None:
        groups = mdata.obs[group_by].unique().tolist()

    for i, cat in enumerate(groups):
        log.info(f"create metacells for group '{cat}' ... ({i+1}/{len(groups)})")
        gc.collect()

        md_sub = mdata.copy()
        md_sub.obsp = None
        mu.pp.filter_obs(md_sub['rna'], group_by, lambda x: x == cat)
        mu.pp.filter_obs(md_sub['atac'], group_by, lambda x: x == cat)
        md_sub.update()
        md_sub = prep_rna(md_sub)
        md_sub = prep_atac(md_sub)
        md_sub = prep_wnn(md_sub)
        # Get metacell groups using WNN graph
        meta_res = get_metacells(md_sub, **kwargs)
        # Add prefix to group numbers
        meta_res['mc_group'] = f"{cat}_" + meta_res['mc_group'].astype(str)
        dfs.append(meta_res)

    # Combine all dataframes
    df = pd.concat(dfs, ignore_index=True)
    return df

