#!/bin/python3
#SBATCH -A kas_dev
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p standard
#SBATCH -t 50:00:00
#SBATCH -o /sfs/lustre/bahamut/scratch/jws6pq/CMfiles/filesort.out
#SBATCH -e /sfs/lustre/bahamut/scratch/jws6pq/CMfiles/filesort.err
#Run program
import os
import sys
sys.path.append('/sfs/lustre/bahamut/scratch/jws6pq/CMfiles')
import pickleopener
FolderDir='/sfs/lustre/bahamut/scratch/jws6pq/Notebook/Finished/'
for file in os.listdir(FolderDir):
    if file[0]=='3':
        pickleopener.GeneratePDBandPlddtfromAlphafold(FolderDir,file,'/sfs/lustre/bahamut/scratch/jws6pq/Notebook/Overall/'+file+'.pdb','/sfs/lustre/bahamut/scratch/jws6pq/Notebook/Overall/'+file+'.plddt')
        pickleopener.GeneratePDBandPlddtfromAlphafold(FolderDir,file,'/sfs/lustre/bahamut/scratch/jws6pq/Notebook/PDB/'+file+'.pdb','/sfs/lustre/bahamut/scratch/jws6pq/Notebook/Plddt/'+file+'.plddt')

