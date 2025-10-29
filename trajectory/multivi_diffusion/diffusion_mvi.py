#!/usr/bin/env python

# conda activate scFates
# nohup python diffusion_mvi.py input_mvi.h5ad out_prefix > logs/diffusion_mvi.log &

import os, sys
import logging as log
import numpy as np
import pandas as pd
import scanpy as sc
import pegasus as pg
import pegasusio as io
import palantir

# Configure logging
log.basicConfig(level=log.INFO)
sc.settings.verbosity = 3
sc.settings.logfile = sys.stdout

# set
new_anno = "anno_wnn_v51"
work_dir = '/work/DevM_analysis/03.trajectory/diffusion_mvi'

# args
indata = sys.argv[1]
output = sys.argv[2]

# read
log.info(f"read input h5ad")
adata = sc.read_h5ad(indata)
log.info(f"anndata info: {adata}")

# palantir
log.info(f"run palantir")
multivi_projections = pd.DataFrame(adata.obsm["X_multiVI"], index=adata.obs_names)
dm_res = palantir.utils.run_diffusion_maps(multivi_projections)
ms_data = palantir.utils.determine_multiscale_space(dm_res)

# neighbor
adata.obsm["X_palantir"] = ms_data.values
sc.pp.neighbors(adata, n_neighbors=30, use_rep="X_palantir")

# draw ForceAtlas2
log.info(f"run draw_graph")
sc.tl.draw_graph(adata)

# PAGA on palantir
log.info(f"run paga")
sc.tl.paga(adata, groups=f"{new_anno}")
sc.pl.paga(adata, color=[new_anno])

# draw ForceAtlas2 with PAGA init_pos
log.info(f"run draw_graph with paga")
sc.tl.draw_graph(adata, init_pos="paga", key_added_ext="fa_paga")

# fle
mdata = io.MultimodalData(adata)
print(mdata)
pg.neighbors(mdata, rep="multiVI")
log.info(f"run diffmap")
pg.diffmap(mdata, rep='multiVI')
pg.fle(mdata, rep="diffmap", target_steps=15000)

# fle with X_palantir
pg.fle(mdata, rep="palantir", out_basis='fle_palantir', target_steps=15000)

# mdata to adata
adata = mdata.to_anndata()

# save
log.info(f"save files...")
adata.write(f'{work_dir}/data/{output}.h5ad')

# draw_graph_fa
dr_df = pd.DataFrame(adata.obsm['X_draw_graph_fa'], index=adata.obs_names, columns=['FA_1', 'FA_2'])
dr_df = dr_df.rename_axis('barcode').reset_index()
dr_df.to_csv(f"{work_dir}/data/{output}_fa.csv", index=False)

# draw_graph_fa_paga
dr_df = pd.DataFrame(adata.obsm['X_draw_graph_fa_paga'], index=adata.obs_names, columns=['FA_PAGA_1', 'FA_PAGA_2'])
dr_df = dr_df.rename_axis('barcode').reset_index()
dr_df.to_csv(f"{work_dir}/data/{output}_fa_paga.csv", index=False)

# fle
dr_df = pd.DataFrame(adata.obsm['X_fle'], index=adata.obs_names, columns=['FLE_1', 'FLE_2'])
dr_df = dr_df.rename_axis('barcode').reset_index()
dr_df.to_csv(f"{work_dir}/data/{output}_fle.csv", index=False)

# fle_palantir
dr_df = pd.DataFrame(adata.obsm['X_fle_palantir'], index=adata.obs_names, columns=['FLE_1', 'FLE_2'])
dr_df = dr_df.rename_axis('barcode').reset_index()
dr_df.to_csv(f"{work_dir}/data/{output}_fle_palantir.csv", index=False)
