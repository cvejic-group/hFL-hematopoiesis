#!/bin/bash
#SBATCH -J HSC_DevDiff
#SBATCH -o logs/HSC_DevDiff_%A_%a.out   # output file %A is job %a is job array index
#SBATCH -a 1-500             # Job Array
#SBATCH -N 1                  # number of nodes
#SBATCH --cpus-per-task=30    # number of threads to execute OpenMP application with
#SBATCH --ntasks-per-node=1   # number of cores per node

source ~/software/anaconda3/bin/activate
conda activate R432

celltype=HSC
work_dir=/work/home/project/20231127_DevM/E2G/lmmPG/$celltype
cd $work_dir

# TIME
input=~/project/20231127_DevM/devm_r432/output/HSC.pg_prep/obj_lst.qs
outdir=$work_dir/PG_DevDiff
if [ ! -d $outdir ]; then mkdir -p $outdir; fi
Rscript ~/project/20231127_DevM/E2G/lmmPG/lmmPG_batch.R -n $SLURM_ARRAY_TASK_ID -i $input -o $outdir --time


