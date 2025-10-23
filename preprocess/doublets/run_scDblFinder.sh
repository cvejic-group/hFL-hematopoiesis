#!/usr/bin/env bash

# nohup bash run_scDblFinder.sh > run_scDblFinder.log &

# 20240123, 40 samples
sed '1d' ~/project/20231127_DevM/sam_info.csv | cut -f 1,4 -d ',' | while read line
do
  sample="$(cut -d ',' -f1 <<< "$line")"
  cellranger="$(cut -d ',' -f2 <<< "$line")"
  #echo $sample
  #echo $cellranger
  Rscript run_scDblFinder.R $sample $cellranger
done

# collect
Rscript collect_doublet_rate.R

