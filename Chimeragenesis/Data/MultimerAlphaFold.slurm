#!/bin/bash
#SBATCH -p gpu          # partition
#SBATCH --gres=gpu:a100:1   # number of GPUs
#SBATCH -N 1           # number of nodes
#SBATCH -n 8            # number of cores
#SBATCH --mem=600000
#SBATCH -t 72:00:00     # time