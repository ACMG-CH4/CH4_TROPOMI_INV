#!/bin/bash

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 0-00:10
#SBATCH -p test
#SBATCH --mem=10000
#SBATCH --mail-type=END

#export OMP_NUM_THREADS=$SLURM_NTASKS
export OMP_NUM_THREADS=1
python Step1_invert.py
exit 0
#EOC
