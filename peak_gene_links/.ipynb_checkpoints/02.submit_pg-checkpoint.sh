#!/bin/bash
#SBATCH --job-name=$1
#SBATCH --output=logs/${1}_%A_%a.out   # output file %A is job %a is job array index
#SBATCH -a 1-500
#SBATCH --partition=small
#SBATCH --account=project_465001709
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=200G
#SBATCH --time=23:59:00
#SBATCH --requeue

source $HOME/software/anaconda3/bin/activate
conda activate R432

celltype=$1
work_dir=~/project/231127_DevM/lmmPG

# ALL
input=$work_dir/$celltype/obj_lst.qs
outdir=$work_dir/$celltype/PG
if [ ! -d $outdir ]; then mkdir -p $outdir; fi
Rscript $work_dir/lmmPG_batch.R -n $SLURM_ARRAY_TASK_ID -i $input -o $outdir --nocorrect

