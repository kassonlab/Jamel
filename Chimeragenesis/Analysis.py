#! python
#
# Code 2023 by Jamel Simpson
import re
from pickle import load as p_load
from json import load as j_load
from os import system, path, strerror, listdir, makedirs
from shutil import copy
from numpy import savetxt, empty, zeros, array, mean
from errno import ENOENT
from Bio import SeqIO, PDB
from collections import defaultdict
from pathlib import Path


def determine_columns_from_container(container):
    data_columns = {}
    # Checks which data columns are wanted by the user by looking for True in the first index of each array in
    # column_names from the analysis_args, Each container can have its own column preferences and every container
    # will have its own columns of data per data column requested
    column_choices = container.analysis_args.column_names
    list_of_chis = container.chimeras
    for data_type, [boolean, title] in column_choices.items():
        if boolean:
            data_columns[title] = tuple(getattr(chimera, data_type) for chimera in list_of_chis)
    return data_columns


def convert_data_dict_to_csv(data_dict, container):
    data_array = empty(((len(container.chimeras) + 1), len(data_dict)), dtype=object)
    for column_count, (column_name, data) in enumerate(data_dict.items()):
        data_array[0, column_count] = column_name
        data_array[1:, column_count] = data
    savetxt(container.analysis_args.analysis_output_csv, data_array, fmt=','.join('%s' for x in data_dict))


def convert_array_to_file(arr, delimiter, file_name):
    savetxt(file_name, arr, fmt=delimiter.join(('%s' for x in arr)))


def get_plddt_dict_from_pdb(pdb_file) -> dict:
    pdb = PDB.PDBParser().get_structure('pdb', pdb_file)[0]
    plddt_dict = defaultdict(tuple)
    chains = tuple(chain for chain in pdb)
    seq_dict = {chain.id.split(':')[-1]: str(chain.seq) for chain in SeqIO.parse(pdb_file, 'pdb-atom')}
    for chain in chains:
        sequence = seq_dict[chain.id]
        plddt = tuple(resi['CA'].bfactor for resi in chain)
        plddt_dict[sequence] += (plddt,)
    plddt_dict = {sequence: tuple(mean(scores) for scores in zip(*plddts)) for sequence, plddts in plddt_dict.items()}
    return plddt_dict


def get_info_from_plddt_file(plddt_file: str):
    with open(plddt_file, 'r') as plddt:
        plddt_text = plddt.read()
        sequence = re.search(r'[a-zA-Z]+', plddt_text).group()
        return sequence, tuple(float(score) for score in re.findall(r'\d+\.\d+', plddt_text))


def create_plddt_file_from_pdb(pdb_file, new_plddt_file):
    with open(new_plddt_file, 'w') as new_plddt:
        for sequence, plddt in get_plddt_dict_from_pdb(pdb_file).items():
            new_plddt.write(
                '>{0}\n{1}'.format(sequence, "\t".join(plddt)))


def skeletonize_alphafold_folder(alphafold_dir, storage_dir):
    rank_file = path.join(alphafold_dir, '../ranking_debug.json')
    makedirs(storage_dir, exist_ok=True)
    new_files = []
    dir_files = [file for file in listdir(alphafold_dir) if file.startswith('ranked')]
    try:
        copy(rank_file, path.join(storage_dir, '../ranking_debug.json'))
        for pdb_file in dir_files:
            new_pdb = path.join(storage_dir, pdb_file)
            new_plddt = path.join(storage_dir, str(Path(pdb_file).stem) + '.plddt')
            new_files.append(new_pdb)
            new_files.append(new_plddt)
            copy(path.join(alphafold_dir, pdb_file), new_pdb)
            create_plddt_file_from_pdb(path.join(alphafold_dir, pdb_file), new_plddt)
        if all(path.exists(file) for file in new_files):
            return True
    except FileNotFoundError:
        print('False')
        return False


def run_Foldx(foldx_file, pdb_file, foldx_command):
    """Call FoldX.  This is optional functionality."""
    pdb_dir = path.dirname(pdb_file)
    foldx_dir = path.dirname(foldx_file)
    system(f'{foldx_command}  -c Stability --pdb {path.basename(pdb_file)} '
           f'--output-dir {foldx_dir} --output-file {path.basename(foldx_file)} '
           f'--pdb-dir {pdb_dir}')


def get_Foldx_results(foldx_file):
    """Read FoldX results."""
    with open(foldx_file, 'r') as foldx_score:
        return foldx_score.read().split()[1]


def generate_alphafold_files(alphafold_folder, output_folder):
    """Creates a text file containing the plddt values of the highest_rank_model extracted from alphafold's result pkl file
    and renames the ranked_0.pdb file and places it in the desired directory."""
    # Checking to see if ranking_debug.json exists. This file is the last to be output by alphafold and is a check that
    # the pkl file you want to extract from exists, as well as to avoid errors
    # ranking_file = path.join(alphafold_folder, '../ranking_debug.json')
    alphafold_folder = Path(alphafold_folder)
    output_folder = Path(output_folder)
    pdb_out = output_folder.joinpath('PDB')
    plddt_out = output_folder.joinpath('Plddt')
    makedirs(pdb_out, exist_ok=True)
    makedirs(plddt_out, exist_ok=True)
    for folder in alphafold_folder.iterdir():
        if (folder.joinpath('ranking_debug.json')).exists():
            # The highest ranked structure is copied with a new name and directory
            old_pdb = folder.joinpath('ranked_0.pdb')
            new_pdb = pdb_out.joinpath(f'{folder.stem}.pdb')
            copy(old_pdb, new_pdb)
            create_plddt_file_from_pdb(new_pdb, plddt_out.joinpath(f'{folder.stem}.plddt'))


def get_sequence_similarity(emboss_file):
    """Returns sequence similarity from an emboss needle file."""
    with open(emboss_file, 'r') as infile:
        infile = infile.read().split('#')
        for line in infile:
            if 'Similarity' in line:
                emboss_score = line.split()[-1].replace('(', '').replace(')', '').replace('%', '')
    return emboss_score


def overall_confidence_from_file(plddt_file):
    """Returns the average confidence score from a protein's plddt file."""
    with open(plddt_file, 'r') as infile:
        plddt = tuple(float(score) for score in infile.readlines())
    average_plddt = sum(plddt) / len(plddt)
    return average_plddt


def overall_confidence(plddt: tuple):
    """Returns the average confidence score from a protein's plddt file."""
    average_plddt = sum(plddt) / len(plddt)
    return average_plddt


def relative_stabilty(plddt_dict: dict[str, tuple], inheritance_dict: dict[str, set], chi_seq):
    raw_stability = 0
    calculations=0
    chimera_label = (set(plddt_dict.keys()) - set(inheritance_dict.keys())).pop()
    for parent, aln_positions in inheritance_dict.items():
        for pos in aln_positions:
            raw_stability += plddt_dict[parent][pos] - plddt_dict[chimera_label][pos]
            calculations+=1
    return raw_stability/calculations
