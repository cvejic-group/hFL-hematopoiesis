'''
conda activate moscot_env
cd /work/DevM_analysis/01.annotation/10.integration_joint_clean_MultiVI/
nohup python 1.2.MultiVI_FL_allPCW.py > 1.2.MultiVI_allPCW.log &
'''

import gzip
import os
import tempfile
from pathlib import Path
import matplotlib.pyplot as plt
import mudata as md
import muon
import numpy as np
#import pooch
import scanpy as sc
import scvi
import seaborn as sns
import torch
import pandas as pd
import string
import anndata as ad

scvi.settings.seed = 0
print("Last run with scvi-tools version:", scvi.__version__)

torch.set_float32_matmul_precision("high")

##---------------------------------------------##
##----------------0. Get ready-----------------##
##---------------------------------------------##

print("Get ready")

adata_paired=sc.read_h5ad("/work/DevM_analysis/01.annotation/10.integration_joint_clean_MultiVI/data/adata_FL_allPCW_concat_RNA-ATAC.h5ad")

##---------------------------------------------##
##--------------1. Prepare object--------------##
##---------------------------------------------##

print("Prepare objects")

adata_mvi=scvi.data.organize_multiome_anndatas(adata_paired)

##---------------------------------------------##
##-----------------2. MultiVI------------------##
##---------------------------------------------##

print("Set up multiVI")

scvi.model.MULTIVI.setup_anndata(
    adata_mvi,
    batch_key="modality", 
    categorical_covariate_keys = ["libraryID", "donorID"],
    layer = "counts",
)

#n_genes=len(rna.var_names)
#n_regions=len(atac.var_names)

mvi=scvi.model.MULTIVI(
    adata_mvi,
    n_genes = 2990,
    n_regions = 50000,
)

mvi.view_anndata_setup()

print("Train multiVI")
mvi.train()

##---------------------------------------------##
##------------------3. Save--------------------##
##---------------------------------------------##

print("Save")
mvi.save("/work/DevM_analysis/01.annotation/10.integration_joint_clean_MultiVI/data/model_multiVI_FL_allPCW",
           overwrite=True, save_anndata=True)
