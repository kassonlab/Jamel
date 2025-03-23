#!python
#
# Script to take AlphaFold output PDBs and set them up for production simulations.
# Code 2022 by Peter Kasson

from os import system, path
from pathlib import Path

def prep(pdb, gmxbin, pdb2gmx):
    shortname = Path(pdb).stem
    system(f'{gmxbin} {pdb2gmx} -f {shortname}.pdb -o {shortname} -p {shortname} -i {shortname}')
    system(f'{gmxbin} editconf -f {shortname}.gro -o {shortname}_box -bt o -c -d 1')
    system(f'{gmxbin} grompp -p {shortname} -c {shortname}_box -f minim.mdp -o {shortname}_em -maxwarn 1')
    system(f'{gmxbin} mdrun -v -deffnm {shortname}_em')
    system(f'{gmxbin} grompp -p {shortname} -c {shortname} -f minim.mdp -o {shortname}_em -maxwarn 1')
    system(f'{gmxbin} solvate -cp {shortname}_em -o {shortname}_sol -p {shortname}')
    system(f'{gmxbin} grompp -p {shortname} -c {shortname}_sol -f minim.mdp -o {shortname}_preion -maxwarn 1')
    system(f'echo 13 | {gmxbin} genion -s {shortname}_preion -o {shortname}_ions -p {shortname} -neutral -conc 0.15')
    system(f'{gmxbin} grompp -p {shortname} -c {shortname}_ions -f minim.mdp -o {shortname}_em -maxwarn 1')
    system(f'{gmxbin} mdrun -v -deffnm {shortname}_em')
    system(f'{gmxbin} grompp -p {shortname} -c {shortname}_em -f startup.mdp -o {shortname}_start -maxwarn 1')
    system(f'{gmxbin} mdrun -v -deffnm {shortname}_start')
    system(f'{gmxbin} grompp -p {shortname} -c {shortname}_start -f prod.mdp -o {shortname}_prod -maxwarn 1')


def create_gmx_setup_slurm(pdb, gmxbin, pdb2gmx, slurm_template, slurm_file_name, output_file, error_file,
                           interpreter='#!/bin/bash'):
    shortname = Path(pdb).stem
    with open(slurm_file_name, 'w') as slurm:
        slurm.write(f'{interpreter}\n#SBATCH -o {output_file}\n#SBATCH -e {error_file}\n')
        with open(slurm_template, "r") as template:
            slurm.write(template.read())
        slurm.write(f'{gmxbin} {pdb2gmx} -f {shortname}.pdb -o {shortname} -p {shortname} -i {shortname}\n')
        slurm.write(f'{gmxbin} editconf -f {shortname}.gro -o {shortname}_box -bt o -c -d 1\n')
        slurm.write(f'{gmxbin} grompp -p {shortname} -c {shortname}_box -f minim.mdp -o {shortname}_em -maxwarn 1\n')
        slurm.write(f'{gmxbin} mdrun -v -deffnm {shortname}_em\n')
        slurm.write(f'{gmxbin} grompp -p {shortname} -c {shortname} -f minim.mdp -o {shortname}_em -maxwarn 1\n')
        slurm.write(f'{gmxbin} solvate -cp {shortname}_em -o {shortname}_sol -p {shortname}\n')
        slurm.write(
            f'{gmxbin} grompp -p {shortname} -c {shortname}_sol -f minim.mdp -o {shortname}_preion -maxwarn 1\n')
        slurm.write(
            f'echo 13 | {gmxbin} genion -s {shortname}_preion -o {shortname}_ions -p {shortname} -neutral -conc 0.15\n')
        slurm.write(f'{gmxbin} grompp -p {shortname} -c {shortname}_ions -f minim.mdp -o {shortname}_em -maxwarn 1\n')
        slurm.write(f'{gmxbin} mdrun -v -deffnm {shortname}_em\n')
        slurm.write(
            f'{gmxbin} grompp -p {shortname} -c {shortname}_em -f startup.mdp -o {shortname}_start -maxwarn 1\n')
        slurm.write(f'{gmxbin} mdrun -v -deffnm {shortname}_start\n')
        slurm.write(f'{gmxbin} grompp -p {shortname} -c {shortname}_start -f prod.mdp -o {shortname}_prod -maxwarn 1')


def create_gmx_prod_slurm(pdb, gmxbin, slurm_template, slurm_file_name, output_file, error_file,
                          interpreter='#!/bin/bash'):
    shortname = Path(pdb).stem
    with open(slurm_file_name, 'w') as slurm:
        slurm.write(f'{interpreter}\n#SBATCH -o {output_file}\n#SBATCH -e {error_file}\n')
        with open(slurm_template, "r") as template:
            slurm.write(template.read())
        slurm.write(f'{gmxbin} mdrun -cpi {shortname}_prod.cpt -v -deffnm {shortname}_prod\n')


# TODO srun and ntomp for mdrun
# TODO sbatch roduction run from inside setup





