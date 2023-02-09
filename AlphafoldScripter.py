from sys import argv
from os import system
from json import load

fastas_to_run = argv[1]
argument_json = argv[2]
with open(fastas_to_run, 'r') as fasta_list:
    fasta_list = fasta_list.readlines()
    fasta_list = [x.split()[0] for x in fasta_list]
with open(argument_json, 'rb') as jfile:
    argument_dict = load(jfile)["arguments"]

proteins_per_slurm = argument_dict['proteins_per_slurm']
file_index = 0
slurm_index = 0
template_slurm = argument_dict['template_slurm'][1]
alphafold_shell_script = argument_dict['alphafold_shell_script'][1]
output_directory = argument_dict['output_directory']
naming_convention = argument_dict['slurm_naming']
character_to_replace = argument_dict['character_to_replace']
while file_index in range(len(fasta_list)):
    current_slurm = naming_convention.replace(character_to_replace, str(slurm_index))
    if argument_dict['template_slurm'][0] == '':
        system(f'cp {template_slurm} {current_slurm}')
        with open(current_slurm, 'a') as slurm_file:
            slurm_file.write(f'\n#SBATCH -o {argument_dict["slurm_output"].replace(character_to_replace, str(slurm_index))}\n'
                             f'#SBATCH -e {argument_dict["slurm_error"].replace(character_to_replace, str(slurm_index))}\n#Run program\n')
            proteins_to_run = ','.join(fasta_list[file_index:file_index + proteins_per_slurm])
            slurm_file.write(f'{alphafold_shell_script} {proteins_to_run} {output_directory}')
    if argument_dict['alphafold_shell_script'][0] == '':
        system(f'sbatch {current_slurm}')
    file_index += proteins_per_slurm
    slurm_index += 1
