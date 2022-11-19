#!/bin/python3
#SBATCH -A kas_dev
#SBATCH -N 1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=25000
#SBATCH -p gpu
#SBATCH --gres=gpu:a100:1
#SBATCH -t 50:00:00
#SBATCH -o /scratch/jws6pq/BridCMfiles/filesort.out
#SBATCH -e /scratch/jws6pq/BridCMfiles/filesort.err
#Run program
from os import listdir
import sys
from concurrent.futures import ProcessPoolExecutor
sys.path.append('/scratch/jws6pq/BridCMfiles')
from pickleopener import GeneratePDBorPlddtfromAlphafold
Files=[x for x in listdir('/scratch/jws6pq/Notebook/Finished/')]
FolderDir=['/scratch/jws6pq/Notebook/Finished/' for x in Files]
Pdb=['/scratch/jws6pq/Notebook/PDB/'+file+'.pdb' for file in Files]
Plddt=['/scratch/jws6pq/Notebook/Plddt/'+file+'.plddt' for file in Files]

with ProcessPoolExecutor(max_workers=25) as executor:
	executor.map(GeneratePDBorPlddtfromAlphafold,FolderDir,Files,Pdb,Plddt)
	Pdb=['/scratch/jws6pq/Notebook/Overall/'+file+'.pdb' for file in Files]
	Plddt=['/scratch/jws6pq/Notebook/Overall/'+file+'.plddt' for file in Files]
	executor.map(GeneratePDBorPlddtfromAlphafold,FolderDir,Files,Pdb,Plddt)

        
        

