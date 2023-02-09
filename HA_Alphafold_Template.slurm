#!/bin/bash
#SBATCH -p gpu          # partition
#SBATCH --gres=gpu   # number of GPUs
#SBATCH -N 1           # number of nodes
#SBATCH -n 6            # number of cores
#SBATCH --mail-user=jws6pq@virginia.edu
#SBATCH --mail-type=END,BEGIN,FAIL
#SBATCH -t 72:00:00     # time