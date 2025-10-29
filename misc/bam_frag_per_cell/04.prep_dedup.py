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

# out dir
outdir=$WORKDIR/bam_per_cell/$cell
""" %(workdir, ppn, cell)
    file_h.write(str2print)


def dedup(file_h=None):
    str2print = """
# deup
$path2samtools view -@ $PPN -h -F 1024 $outdir/${cell}.sorted.bam | $path2samtools sort -@ $PPN -O bam -o $outdir/${cell}.dedup.bam
$path2samtools index -@ $PPN $outdir/${cell}.dedup.bam
"""
    file_h.write(str2print)


def dedup_pcw(file_h=None, pcw=None):
    str2print = """
pcw=%s
$path2samtools view -@ $PPN -h -F 1024 $outdir/${cell}_PCW${pcw}.sorted.bam | $path2samtools sort -@ $PPN -O bam -o $outdir/${cell}_PCW${pcw}.dedup.bam
$path2samtools index -@ $PPN $outdir/${cell}_PCW${pcw}.dedup.bam
""" %(pcw)
    file_h.write(str2print)


def run():
    """The main function
    """
    workdir = "/work/home/project/20231127_DevM/sinto"
    outdir = f"{workdir}/bam_per_cell"
    PPN = 24
    # cells
    cells = ["HSC", "GP", "Granulocyte",
             "MEMP-t", "MEMP", "MEP", "MEMP-Mast-Ery", "MEMP-Ery", "Early-Ery", "Late-Ery",
             "MEMP-MK", "MK", "MastP-t", "MastP", "Mast",
             "MDP", "Monocyte", "Kupffer", "cDC1", "cDC2", "pDC", "ASDC",
             "LMPP", "LP", "Cycling-LP", "PreProB", "ProB-1", "ProB-2", "Large-PreB", "Small-PreB", "IM-B",
             "NK", "ILCP", "T",
             "Hepatocyte", "Endothelia"]
    for cell in cells:
        bash_file = f"{outdir}/{cell}/{cell}.dedup.sh"
        file_h = open(bash_file, 'w')
        header_info(file_h, workdir, PPN, cell)
        dedup(file_h)
        # per pcw
        for pcw in [str(i) for i in range(5, 19)]:
            if os.path.exists(f"{outdir}/{cell}/{cell}_PCW{pcw}.sorted.bam"):
                dedup_pcw(file_h, pcw)
        file_h.close()


if __name__ == '__main__':
    run()
