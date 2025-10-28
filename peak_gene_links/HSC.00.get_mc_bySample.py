#!/usr/bin/env python

# conda activate muon
# nohup python get_hsc_mc_bySample.py &

import muon as mu
import pandas as pd
from get_metacell_wnn import *

# read
mdata = mu.read(f"/work/DevM_analysis/01.annotation/11.subclustering/HSC/data/FL_wnn_clustered.v00.h5mu")
mdata

# by sample (per donor per lib)
df_mc = get_metacells_by_group(mdata=mdata, group_by='sampleID')

# save
df_mc.to_csv("hsc_mc_bySample.csv", index=False)

