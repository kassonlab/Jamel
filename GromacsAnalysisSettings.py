#!/bin/python
#SBATCH -p gpu          # partition
#SBATCH --gres=gpu  # number of GPUs
#SBATCH -N 1            # number of nodes
#SBATCH --ntasks-per-node=4
#SBATCH --mail-user=jws6pq@virginia.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -t 72:00:00     # time
#SBATCH -A kassonlab
#SBATCH -o /scratch/jws6pq/Gromacs/analysis.out
#SBATCH -e /scratch/jws6pq/Gromacs/analysis.err
#Run program
import sys
from numpy import savetxt, array
from concurrent.futures import ProcessPoolExecutor
sys.path.append("/gpfs/gpfs0/scratch/jws6pq/Gromacs/")
from Chimeragenesis.Scripts import GromacsAnalysis

GromacsAnalysis.call_gromacs_in_rivanna()
protein_list= [x.split()[0] for x in open("/gpfs/gpfs0/scratch/jws6pq/Gromacs/List", 'r').readlines()][8:]
tprfile=[f"/gpfs/gpfs0/scratch/jws6pq/Gromacs/3merSARS2w{x}S1_production_1.tpr" for x in protein_list]
xtcfile=[f"/gpfs/gpfs0/scratch/jws6pq/Gromacs/3merSARS2w{x}S1_production_1.xtc" for x in protein_list]
timestep=['1000' for x in protein_list]
newxtcfile=[f"/gpfs/gpfs0/scratch/jws6pq/Gromacs/3merSARS2w{x}dt{timestep[0]}.xtc" for x in protein_list]
xvgfile=[f"/gpfs/gpfs0/scratch/jws6pq/Gromacs/3merSARS2w{x}dt{timestep[0]}.xvg" for x in protein_list]
averagedxvgfile=[f"/gpfs/gpfs0/scratch/jws6pq/Gromacs/Avg3merSARS2w{x}dt{timestep[0]}.xvg" for x in protein_list]
movie_pdb=[f"/gpfs/gpfs0/scratch/jws6pq/Gromacs/3merSARS2w{x}S1movie.pdb" for x in protein_list]
rmsd_time_steps=[str(time) for time in range(0,1001,100)]
with ProcessPoolExecutor() as exe:
    exe.map(GromacsAnalysis.truncate_xtc_file, tprfile, xtcfile, newxtcfile, timestep)
with ProcessPoolExecutor() as exe:
    exe.map(GromacsAnalysis.create_rmsf_file, timestep, tprfile, newxtcfile, xvgfile)
with ProcessPoolExecutor() as exe:
    exe.map(GromacsAnalysis.averaging_multimer_rmsf, xvgfile, averagedxvgfile)
with ProcessPoolExecutor() as exe:
    exe.map(GromacsAnalysis.create_trajectory_movie_pdb, tprfile, newxtcfile, movie_pdb)
all_rmsd_changes=[]
for protein_index, protein in enumerate(protein_list):
    rmsd_change=[protein]
    for time in rmsd_time_steps:
        GromacsAnalysis.create_pdb_from_trajectory(tprfile[protein_index], newxtcfile[protein_index],
                                                   f'/gpfs/gpfs0/scratch/jws6pq/Gromacs/SARSw{protein}S1frame{time}.pdb', time)
        rmsd_change.append(
            GromacsAnalysis.pymol_rmsd(f'/gpfs/gpfs0/scratch/jws6pq/Gromacs/SARSw{protein}S1frame{time}.pdb',
                                                      f'/gpfs/gpfs0/scratch/jws6pq/Gromacs/SARSw{protein}S1frame0.pdb'))
    all_rmsd_changes.append(rmsd_change)

savetxt('/gpfs/gpfs0/scratch/jws6pq/Gromacs/122022ChangeinRMSD.tsv', array(all_rmsd_changes, dtype=object),
        fmt='%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s')




