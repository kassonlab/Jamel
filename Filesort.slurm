#!/bin/python3
#SBATCH -A kas_dev
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=25000
#SBATCH -p gpu
#SBATCH --gres=gpu
#SBATCH -t 50:00:00
#SBATCH -o /scratch/jws6pq/BridCMfiles/filesort.out
#SBATCH -e /scratch/jws6pq/BridCMfiles/filesort.err
#Run program
from os import listdir
import sys
from concurrent.futures import ProcessPoolExecutor
sys.path.append('/gpfs/gpfs0/scratch/jws6pq/Storage/')
from pickleopener import generate_alphafold_files,limited_alphafold_transfer
Files=[x for x in listdir('/scratch/jws6pq/Notebook/Finished/')]
FolderDir=['/scratch/jws6pq/Notebook/Finished/' for x in Files]
Pdb=[f'/scratch/jws6pq/Notebook/PDB/{file}.pdb' for file in Files]
Plddt=[f'/scratch/jws6pq/Notebook/Plddt/{file}.plddt' for file in Files]

with ProcessPoolExecutor(max_workers=25) as executor:
	executor.map(generate_alphafold_files, FolderDir, Files, Pdb, Plddt)
	Pdb=[f'/scratch/jws6pq/Notebook/Overall/{file}.pdb' for file in Files]
	Plddt=[f'/scratch/jws6pq/Notebook/Overall/{file}.plddt' for file in Files]
	executor.map(generate_alphafold_files, FolderDir, Files, Pdb, Plddt)
	executor.map(limited_alphafold_transfer, FolderDir, Files, ['/gpfs/gpfs0/scratch/jws6pq/Storage/' for file in Files],
				 ['/gpfs/gpfs0/scratch/jws6pq/Storage/' for file in Files])

        

