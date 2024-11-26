#!python
#
# Script to take AlphaFold output PDBs and set them up for production simulations.
# Code 2022 by Peter Kasson

from glob import glob
from os import system, path
from pathlib import Path

SLURM_TEMPLATE = 'MultimerAlphaFold.slurm'
SH_TEMPLATE = 'MultimerAlphaFold.sh'


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

def create_alphafold_slurm(iter_of_fastas, slurm_filename, output_directory,subunit_count):
    system(f'cp {SLURM_TEMPLATE} {slurm_filename}')
    subunits = 'monomer' if subunit_count==1 else 'multimer'
    databases = 'full_dbs' if subunit_count == 1 else 'reduced_dbs'
    with open(slurm_filename, 'a') as slurm_file:
        output_file=Path(output_directory).joinpath(Path(slurm_filename).stem+'.out')
        slurm_file.write(
            f'\n#SBATCH -o {Path(output_directory).joinpath(Path(slurm_filename).stem+".out")}\n'
            f'#SBATCH -e {Path(output_directory).joinpath(Path(slurm_filename).stem+".err")}\n#Run program\n')
        proteins_to_run = ','.join(iter_of_fastas)
        slurm_file.write(f'{SH_TEMPLATE} {proteins_to_run} {subunits} {databases} {output_directory}')


def alphafold_submission_for_chimera_container(container, list_of_fastas):
    # fasta list is a separate input for control
    fasta_to_run = ()
    placeholder = container.naming_args.PLACEHOLDER
    output_directory = container.naming_args.alphafold_outputs_dir
    submission_toggles = container.submission_args.submission_toggles
    proteins_per_slurm = container.submission_args.proteins_per_slurm
    naming_convention = container.submission_args.slurm_naming
    # Loops through all fastas created and checks if they are complete by looking for their ranking_debug file
    if submission_toggles['stragglers_or_custom_or_all'] == 'stragglers':
        for fasta in list_of_fastas:
            if not path.exists(output_directory + Path(fasta).stem + '/ranking_debug.json'):
                fasta_to_run += (fasta,)
        # Puts all fastas in a line separated file specified by custom_list_to_run
        if submission_toggles['create_file_of_stragglers']:
            with open(container.submission_args.custom_list_to_run, 'w') as run_list:
                run_list.write('\n'.join(fasta for fasta in fasta_to_run))
        # if all of them are complete and fasta_to_run is empty then all slurm actions are toggled off
        if not fasta_to_run:
            submission_toggles['create_slurms'] = submission_toggles['sbatch slurms'] = False
            print('All Done')
    elif submission_toggles['stragglers_or_custom_or_all'] == 'custom':
        with open(container.submission_args.custom_list_to_run, 'r') as run_list:
            run_list = run_list.readlines()
            fasta_to_run = [x.split()[0] for x in run_list]
    elif submission_toggles['stragglers_or_custom_or_all'] == 'all':
        fasta_to_run = list_of_fastas
    for slurm_index, file_index in enumerate(range(0, len(fasta_to_run), proteins_per_slurm)):
        current_slurm = naming_convention.replace(placeholder, str(slurm_index))
        if submission_toggles['create_slurms']:
            create_alphafold_slurm(fasta_to_run[file_index:file_index + proteins_per_slurm], current_slurm,
                                   output_directory)
        if submission_toggles['sbatch_slurms']:
                system(f'sbatch {current_slurm}')
