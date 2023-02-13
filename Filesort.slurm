
from os import listdir,path
from json import load
from sys import  argv
from concurrent.futures import ProcessPoolExecutor
from pickleopener import generate_alphafold_files, limited_alphafold_transfer
from pathlib import Path
argument_json = argv[1]
with open(argument_json, 'rb') as jfile:
    argument_dict = load(jfile)["arguments"]

replacing=argument_dict['character_to_replace']
Files = [x for x in listdir(argument_dict['alphafold_directory'])]
alphafold_directory = [argument_dict['alphafold_directory'] for x in Files]
pdb_naming = [argument_dict['pdb_naming'].replace(replacing,file) for file in Files]
plddt_naming = [argument_dict['plddt_naming'].replace(replacing,file) for file in Files]

with ProcessPoolExecutor(max_workers=25) as executor:
    executor.map(generate_alphafold_files, alphafold_directory, plddt_naming,pdb_naming)
    # executor.map(limited_alphafold_transfer, alphafold_directory, Files,
    #              ['/gpfs/gpfs0/scratch/jws6pq/Storage/' for file in Files],
    #              ['/gpfs/gpfs0/scratch/jws6pq/Storage/' for file in Files])
if argument_dict['fasta_list'][0]=='':
    with open(argument_dict['fasta_list'][1], 'r') as list_to_check:
        list_to_check=list_to_check.readlines()
        for fasta in list_to_check:
            if not path.exists(argument_dict['plddt_naming'].replace(replacing, Path(fasta).stem)):
                with open(argument_dict['list_of_stragglers'], 'w') as run_list:
                    run_list.write(f'{fasta}\n')


