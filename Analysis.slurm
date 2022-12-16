#!/bin/bash
#SBATCH -p gpu          # partition
#SBATCH --gres=gpu:a100:1   # number of GPUs
#SBATCH -N 1           # number of nodes
#SBATCH --ntasks-per-node=4       # number of cores
#SBATCH --mail-user=jws6pq@virginia.edu
#SBATCH --mail-type=END,BEGIN,FAIL
#SBATCH -t 72:00:00     # time
#SBATCH -A kas_dev
#Run program
cd /gpfs/gpfs0/scratch/jws6pq/Notebook/Plddt
/gpfs/gpfs0/scratch/jws6pq/pymol/pymol -c AnalysisSetings.py