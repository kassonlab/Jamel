#!python
#
# Script to take AlphaFold output PDBs and set them up for production simulations.
# Code 2022 by Peter Kasson

from glob import glob
from os import system
from pathlib import Path
def prep(pdb,gmxbin,pdb2gmx):
    shortname = Path(pdb).stem
    system(f'{gmxbin} {pdb2gmx} -f {shortname} -o {shortname} -p {shortname} -i {shortname}')
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
# def create_setup_slurm(pdb,gmxbin,pdb2gmx,slurm_template,slurm_file_name):

    # with open()
if __name__ == '__main__':
  for pdbfile in glob.glob('3mer*pdb'):
      prep(pdbfile)
