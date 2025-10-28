#!/usr/bin/env bash

# prep peaks
for cell in NK
do
  bash 01.prep_gene_annot.sh $cell
done

# magma
# Rheumatoid arthritis
bash 02.run_magma.sh NK ra_GCST90132223 97173
bash 02.run_magma.sh NK ra_finn 331429

# dermatitiseczema_finn
bash 02.run_magma.sh NK dermatitiseczema_finn 500348


