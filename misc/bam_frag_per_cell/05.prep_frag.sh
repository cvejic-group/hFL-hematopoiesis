#!/usr/bin/env bash

# prep frag files

cells=$1

for cell in "${cells[@]}"
do

echo $cell
bash frag_per_cell/$cell/$cell.prep_frag.sh

done


