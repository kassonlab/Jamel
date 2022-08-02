#!/usr/bin/python3
#SBATCH -p gpu          # partition
#SBATCH --gres=gpu:v100:1    # number of GPUs
#SBATCH -N 1            # number of nodes
#SBATCH -c 8            # number of cores
#SBATCH -t 72:00:00     # time

import argparse
import os

def run_alpha(seqfile, outputdir):
  """Wrapper to run AlphaFold on one sequence.
  Args:
    seqfile: full path to FASTA file
    outputdir:  full path for output directory."""
  os.system('/scratch/jws6pq/CMfiles/alpha_run.sh %s %s' % (seqfile, outputdir))

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('-s', help='comma-separated list of sequence files',
                     type=str)
  parser.add_argument('-o', help='output directory', type=str)
  args = parser.parse_args()
  seqlist = args.s.split(',')
  for seq in seqlist:
    run_alpha(os.path.abspath(seq), os.path.abspath(args.o))
