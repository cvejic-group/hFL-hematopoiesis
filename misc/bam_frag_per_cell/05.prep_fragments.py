#!/usr/bin/env python

import os
import sys
import numpy as np
import pandas as pd


def header_info(file_h=None, workdir=None, ppn=None, cell=None):
    """header
    """
    str2print = """#!/usr/bin/env bash

# WORK DIR
WORKDIR=%s

# PPN
PPN=%s

# cell
cell=%s

# tools
TOOLDIR=$HOME/software
path2samtools=$TOOLDIR/samtools/1.18/samtools

# sinto
source $TOOLDIR/anaconda3/bin/activate
conda activate sinto

# input dir
inputdir=$WORKDIR/bam_per_cell/$cell

# out dir
outdir=$WORKDIR/frag_per_cell/$cell
""" %(workdir, ppn, cell)
    file_h.write(str2print)


def prep_frag(file_h=None):
    str2print = """
# make frag
sinto fragments --collapse_within -p $PPN -b $inputdir/${cell}.sorted.bam -f $outdir/${cell}.fragments.tmp.tsv
sort -k 1,1 -k2,2n $outdir/${cell}.fragments.tmp.tsv > $outdir/${cell}.fragments.tsv
bgzip -@ $PPN $outdir/${cell}.fragments.tsv
tabix -p bed $outdir/${cell}.fragments.tsv.gz
rm $outdir/${cell}.fragments.tmp.tsv
"""
    file_h.write(str2print)


def prep_frag_pcw(file_h=None, pcw=None):
    str2print = """
pcw=%s
sinto fragments --collapse_within -p $PPN -b $inputdir/${cell}_PCW${pcw}.sorted.bam -f $outdir/${cell}_PCW${pcw}.fragments.tmp.tsv
sort -k 1,1 -k2,2n $outdir/${cell}_PCW${pcw}.fragments.tmp.tsv > $outdir/${cell}_PCW${pcw}.fragments.tsv
bgzip -@ $PPN $outdir/${cell}_PCW${pcw}.fragments.tsv
tabix -p bed $outdir/${cell}_PCW${pcw}.fragments.tsv.gz
rm $outdir/${cell}_PCW${pcw}.fragments.tmp.tsv
""" %(pcw)
    file_h.write(str2print)


def run():
    """The main function
    """
    workdir = "/work/DevM_analysis/data/sinto"
    datadir = f"{workdir}/bam_per_cell"
    outdir = f"{workdir}/frag_per_cell"
    PPN = 22
    # cells
    cells = ["HSC", "GP", "Granulocyte",
             "MEMP-t", "MEMP", "MEP", "MEMP-Mast-Ery", "MEMP-Ery", "Early-Ery", "Late-Ery",
             "MEMP-MK", "MK", "MastP-t", "MastP", "Mast",
             "MDP", "Monocyte", "Kupffer", "cDC1", "cDC2", "pDC", "ASDC",
             "LMPP", "LP", "Cycling-LP", "PreProB", "ProB-1", "ProB-2", "Large-PreB", "Small-PreB", "IM-B",
             "NK", "ILCP", "T",
             "Hepatocyte", "Endothelia"]
    for cell in cells:
        if not os.path.exists(f"{outdir}/{cell}"):
            os.makedirs(f"{outdir}/{cell}", exist_ok=True)
        bash_file = f"{outdir}/{cell}/{cell}.prep_frag.sh"
        file_h = open(bash_file, 'w')
        header_info(file_h, workdir, PPN, cell)
        prep_frag(file_h)
        # per pcw
        for pcw in [str(i) for i in range(5, 19)]:
            if os.path.exists(f"{datadir}/{cell}/{cell}_PCW{pcw}.sorted.bam"):
                prep_frag_pcw(file_h, pcw)
        file_h.close()


if __name__ == '__main__':
    run()

