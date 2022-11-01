#!/bin/bash
#SBATCH -p gpu          # partition
#SBATCH --gres=gpu:a100:1    # number of GPUs
#SBATCH -N 1            # number of nodes
#SBATCH -c 8            # number of cores
#SBATCH --mail-user=jws6pq@virginia.edu
#SBATCH --mail-type=END,FAIL
#SBATCH -t 72:00:00     # time
#SBATCH -A kas_dev
#SBATCH -o /scratch/jws6pq/Gromacs/1gromacs.out
#SBATCH -e /scratch/jws6pq/Gromacs/1gromacs.err
#Run program
module load gcc/9.2.0
module load cuda/11.0.228
module load openmpi/3.1.6
module load gcccuda/9.2.0_11.0.228
module load goolfc/9.2.0_3.1.6_11.0.228
module load gromacs/2021.2
cd /scratch/jws6pq/Gromacs/
gmx_mpi pdb2gmx -f 3merSARS2wBAT2006S1.pdb -ff charmm36-jul2021 -o 3merSARS2wBAT2006S1.gro -p 3merSARS2wBAT2006S1.top -i 3merSARS2wBAT2006S1.itp -water spce -ignh
gmx_mpi editconf -f 3merSARS2wBAT2006S1.gro -o 3merSARS2wBAT2006S1_newbox.gro -c -d 1.0 -bt cubic
gmx_mpi solvate -cp 3merSARS2wBAT2006S1_newbox.gro -cs spc216.gro -o 3merSARS2wBAT2006S1_solv.gro -p 3merSARS2wBAT2006S1.top
gmx_mpi grompp -f ions.mdp -c 3merSARS2wBAT2006S1_solv.gro -p 3merSARS2wBAT2006S1.top -o 3merSARS2wBAT2006S1_ions.tpr -maxwarn 1
echo SOL | gmx_mpi genion -s 3merSARS2wBAT2006S1_ions.tpr -o 3merSARS2wBAT2006S1_ions.gro -p 3merSARS2wBAT2006S1.top -pname SOD -nname CLA -neutral
gmx_mpi grompp -f minim.mdp -c 3merSARS2wBAT2006S1_ions.gro -p 3merSARS2wBAT2006S1.top -o 3merSARS2wBAT2006S1em.tpr
gmx_mpi mdrun -v -deffnm 3merSARS2wBAT2006S1em
gmx_mpi grompp -f nvt.mdp -c 3merSARS2wBAT2006S1em.gro -r 3merSARS2wBAT2006S1em.gro -p 3merSARS2wBAT2006S1.top -o 3merSARS2wBAT2006S1nvt.tpr
gmx_mpi mdrun -deffnm 3merSARS2wBAT2006S1nvt
gmx_mpi grompp -f npt.mdp -c 3merSARS2wBAT2006S1nvt.gro -r 3merSARS2wBAT2006S1nvt.gro -t 3merSARS2wBAT2006S1nvt.cpt -p 3merSARS2wBAT2006S1.top -o 3merSARS2wBAT2006S1npt.tpr
gmx_mpi mdrun -deffnm 3merSARS2wBAT2006S1npt
gmx_mpi grompp -f charmmmd.mdp -c 3merSARS2wBAT2006S1npt.gro -t 3merSARS2wBAT2006S1npt.cpt -p 3merSARS2wBAT2006S1.top -o 3merSARS2wBAT2006S1_production_1.tpr
gmx_mpi mdrun -ntomp 64 -cpi 3merSARS2wBAT2006S1_production_1.cpt -deffnm 3merSARS2wBAT2006S1_production_1