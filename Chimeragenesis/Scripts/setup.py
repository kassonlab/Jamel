#
# Script to take AlphaFold output PDBs and set them up for production simulations.
# Code 2022 by Peter Kasson
import os
import subprocess
import sys
from os import system, path
from pathlib import Path
from ChimeraClasses import HomomerChimeraArgs
from ESM import SequenceDataframe
def submit_multifasta_alphafold(arg_json=''):
    arg_json = HomomerChimeraArgs(sys.argv[1],False) if not arg_json else HomomerChimeraArgs(arg_json,False)
    shell_script=arg_json.submission_args.alphafold_shell_script
    output_direc = arg_json.fasta_args.output_directory
    proteins_per_slurm = arg_json.submission_args.proteins_per_slurm
    aln_df=SequenceDataframe(arg_json.fasta_args.collective_fasta)
    fasta_path=Path(output_direc).joinpath('Fasta')
    alphafold_path=Path(output_direc).joinpath('Alphafold')
    if not Path.exists(fasta_path):
        os.makedirs(fasta_path)
    if not Path.exists(alphafold_path):
        os.makedirs(alphafold_path)
    aln_df.make_individual_fasta(fasta_path,arg_json.fasta_args.number_of_subunits)
    fasta_to_run=[str(file) for file in fasta_path.iterdir() if file.is_file()]
    for slurm_index, file_index in enumerate(range(0, len(fasta_to_run), proteins_per_slurm)):
        command_list=['sbatch', '-J', str(slurm_index) + 'Alphafold' ,
         '-e', Path(output_direc).joinpath(
            'Alphafold' + str(slurm_index)).with_suffix(".err"),
         shell_script,
         ",".join(fasta_to_run[file_index:file_index + proteins_per_slurm]),
         alphafold_path]
        result = subprocess.run(command_list, capture_output=True)
        print(result.stdout,result.stderr)
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

if __name__=='__main__':
    if sys.argv[1]:
        submit_multifasta_alphafold()






