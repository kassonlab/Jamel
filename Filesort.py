
from os import listdir,path
from json import load
from sys import  argv
from concurrent.futures import ProcessPoolExecutor
from pickleopener import generate_alphafold_files, limited_alphafold_transfer
from pathlib import Path
argument_json = argv[1]
with open(argument_json, 'rb') as jfile:
    argument_dict = load(jfile)["arguments"]

with open(argument_dict['list_of_stragglers'], 'w') as run_list:
    with open(argument_dict['fasta_list'], 'r') as list_to_check:
        list_to_check=list_to_check.readlines()
        for fasta in list_to_check:
            print(path.exists(argument_dict['alphafold_directory'] + Path(fasta).stem + '/ranking_debug.json'))
            if not path.exists(argument_dict['alphafold_directory']+Path(fasta).stem+'/ranking_debug.json'):
                run_list.write(f'{fasta}')


