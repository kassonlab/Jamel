#!/bin/python3
#SBATCH -p gpu          # partition
#SBATCH --gres=gpu:a100:1
#SBATCH -N 1           # number of nodes
#SBATCH --ntasks-per-node=7
#SBATCH
#SBATCH --mail-type=ALL
#SBATCH -t 72:00:00     # time
#SBATCH -A
#SBATCH -o /scratch/jws6pq/BridCMfiles/transferpy.out
#SBATCH -e /scratch/jws6pq/BridCMfiles/transferpy.err
#Run program
import os
import argparse
import sys
sys.path.append('/gpfs/gpfs0/scratch/jws6pq/Storage/')
from Analysis import limited_alphafold_transfer
parser = argparse.ArgumentParser(description='')
parser.add_argument('-t', '--transfer', type=str, required=True, help='')
args = parser.parse_args()

# Filestotransfer=['zip -m -r /gpfs/gpfs0/scratch/jws6pq/December2022/'+x+'.zip /gpfs/gpfs0/scratch/jws6pq/Notebook/Finished/'+x for x in os.listdir('/gpfs/gpfs0/scratch/jws6pq/Notebook/Finished/')]
# with ProcessPoolExecutor() as exe:
#  exe.map(os.system,Filestotransfer)


files = [file.replace('.zip', '') for file in os.listdir('/gpfs/gpfs0/scratch/jws6pq/December2022/') if
         file.endswith('zip') if file.startswith('3')]
Filestotransfer = ['unzip -j  /gpfs/gpfs0/scratch/jws6pq/from_old_scratch/OldAlphaOutput/{0}.zip -d /gpfs/gpfs0/scratch/jws6pq/Notebook/Finished/{0}'.format(x) for x in files]

Files = [x for x in os.listdir('/scratch/jws6pq/Notebook/Finished/')]
FolderDir = ['/scratch/jws6pq/Notebook/Finished/' for x in Files]
Pdb = [f'/scratch/jws6pq/Notebook/PDB/{x}.pdb' for x in Files]
Plddt = [f'/scratch/jws6pq/Notebook/Plddt/{x}.plddt' for x in Files]
with ProcessPoolExecutor(max_workers=4) as executor:
    executor.map(limited_alphafold_transfer, FolderDir, Files,
                 ['/gpfs/gpfs0/scratch/jws6pq/Storage/' for file in Files],
                 ['/gpfs/gpfs0/scratch/jws6pq/Storage/' for file in Files])
    executor.map(generate_alphafold_files, FolderDir, Files, Pdb, Plddt)
    Pdb = [f'/scratch/jws6pq/Notebook/Overall/{x}.pdb' for x in Files]
    Plddt = [f'/scratch/jws6pq/Notebook/Overall/{x}.plddt' for x in Files]
    executor.map(generate_alphafold_files, FolderDir, Files, Pdb, Plddt)
for file in Files:
    if os.path.exists(f'/gpfs/gpfs0/scratch/jws6pq/Storage/{file}_4_plddt.npy') and os.path.exists(
                    f'/gpfs/gpfs0/scratch/jws6pq/Storage/{file}_2_plddt.npy') and os.path.exists(
                    f'/gpfs/gpfs0/scratch/jws6pq/Storage/{file}_0_plddt.npy'):
        os.system(f'rm /gpfs/gpfs0/scratch/jws6pq/December2022/{file}.zip')
        os.system(f'rm -r /gpfs/gpfs0/scratch/jws6pq/Notebook/Finished/{file}')
