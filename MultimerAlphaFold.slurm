#!/bin/bash
#SBATCH -p gpu          # partition
#SBATCH --gres=gpu:a100:1   # number of GPUs
#SBATCH -N 1           # number of nodes
#SBATCH -n 10            # number of cores
#SBATCH --mail-user=jws6pq@virginia.edu
#SBATCH --mail-type=END,BEGIN,FAIL
#SBATCH --mem=600000
#SBATCH -t 72:00:00     # time
#SBATCH -A kas_dev