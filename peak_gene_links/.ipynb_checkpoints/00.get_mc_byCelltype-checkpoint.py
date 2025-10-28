#!/usr/bin/env python

# conda activate muon
# nohup python get_mc_byCelltype.py > logs/get_mc_byCelltype.log &

import muon as mu
import pandas as pd
from get_metacell_wnn import *

# read
mdata = mu.read(f"/work/DevM_analysis/01.annotation/11.subclustering/blood/data/FL_wnn_clustered.v00.h5mu")
mdata

# by celltype
df_mc = get_metacells_by_group(mdata=mdata, group_by='anno_wnn_v51')

# save
df_mc.to_csv("blood_celltype_mc.csv", index=False)

