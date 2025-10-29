#!/usr/bin/env bash

cells=$1

for cell in "${cells[@]}"
do

echo $cell
bash bam_per_cell/$cell/$cell.dedup.sh

done


