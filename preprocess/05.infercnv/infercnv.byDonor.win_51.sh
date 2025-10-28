#!/usr/bin/env bash

# nohup bash analysis/infercnv.byDonor.win_51.sh > analysis/infercnv.byDonor.win_51.log &

# ulimit
ulimit -s unlimited

for patient in "${patients[@]}"
do
  echo $patient
  Rscript analysis/infercnv.byDonor.win_51.R $patient
done

# collect results
out_dir="infercnv.byDonor.win_51"
mkdir output/$out_dir/infercnv.byDonor.win_51
for patient in "${patients[@]}"
do
  echo $patient
  cp output/$out_dir/$patient/infercnv.pdf output/$out_dir/infercnv.byDonor.win_51/$patient.pdf
done

