#!/usr/bin/env bash

cells=$1

for cell in "${cells[@]}"
do

echo $cell
bash bam_per_cell/$cell/$cell.sh > bam_per_cell/$cell/$cell.log 2>&1

done


