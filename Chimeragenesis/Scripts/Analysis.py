#! python
#
# Code 2023 by Jamel Simpson
import re
from os import system, path, listdir, makedirs
from shutil import copy
import numpy as np
from numpy import savetxt, empty, mean
from Bio import SeqIO, PDB
from collections import defaultdict
from pathlib import Path
from pandas import DataFrame

from matplotlib import pyplot as plt

def scatterplot_for_df(df:DataFrame, x, y, attr_mods:dict=None, scatter_args:dict=None):
    if scatter_args is None:
        scatter_args = {}
    df.plot.scatter(x=x, y=y,**scatter_args)
    if attr_mods:
        for attr,mod in attr_mods.items():
            plt.attr(mod)

    plt.show()
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


def convert_data_dict_to_csv(data: dict, container):
    data_array = empty(((len(container.chimeras) + 1), len(data)), dtype=object)
    for column_count, (column_name, data) in enumerate(data.items()):
        data_array[0, column_count] = column_name
        data_array[1:, column_count] = data
    savetxt(container.analysis_args.analysis_output_csv, data_array, fmt=','.join('%s' for _ in data))


def convert_array_to_file(arr, delimiter, file_name):
    savetxt(file_name, arr, fmt=delimiter.join(('%s' for _ in arr)))


def get_plddt_dict_from_pdb(pdb_file) -> dict[str, tuple]:
    pdb = PDB.PDBParser().get_structure('pdb', pdb_file)[0]
    plddt_dict = defaultdict(tuple)
    seq_dict = {chain.id.split(':')[-1]: str(chain.seq) for chain in SeqIO.parse(pdb_file, 'pdb-atom')}
    for chain in pdb:
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
                f'>{0}\n{1}'.format(sequence, "\t".join(str(plddt))))


def skeletonize_alphafold_folder(alphafold_dir, storage_dir):
    """Takes an alphafold output folder and copies the ranking_debug,the top ranked pdb and the corresponding plddt scores  into another directory."""
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
    pdb_dir = path.dirname(pdb_file)
    foldx_dir = path.dirname(foldx_file)
    system(f'{foldx_command}  -c Stability --pdb {path.basename(pdb_file)} '
           f'--output-dir {foldx_dir} --output-file {path.basename(foldx_file)} '
           f'--pdb-dir {pdb_dir}')


def get_Foldx_results(foldx_file):
    with open(foldx_file, 'r') as foldx_score:
        return foldx_score.read().split()[1]


def generate_alphafold_files(alphafold_folder, output_folder, make_plddt=False):
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
            if make_plddt:
                create_plddt_file_from_pdb(new_pdb, plddt_out.joinpath(f'{folder.stem}.plddt'))


def get_emboss_sequence_identity(emboss_file):
    """Returns sequence similarity from an emboss needle file."""
    with open(emboss_file, 'r') as infile:
        infile = infile.read()
    emboss_score = re.search(r'# Identity:.*\(([\d.]+)%', infile).group(1)
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


def relative_stabilty(plddt_dict: dict[str, tuple], inheritance_dict: dict[str, dict[tuple, tuple]]):
    raw_stability = []
    chimera_label = (set(plddt_dict.keys()) - set(inheritance_dict.keys())).pop()
    for parent, boundaries in inheritance_dict.items():
        for parent_bound, chimera_bound in boundaries.items():
            if not all(bound == 0 for bound in parent_bound):
                for chi, native in zip(plddt_dict[chimera_label][slice(*chimera_bound)],
                                       plddt_dict[parent][slice(*parent_bound)]):
                    raw_stability.append((chi - native) / native * 100)
    return np.mean(raw_stability)

def relative_stabilty_v2(plddt_dict: dict[str, tuple], inheritance_dict: dict[str, dict[tuple, tuple]]):
    raw_stability = []
    chimera_label = (set(plddt_dict.keys()) - set(inheritance_dict.keys())).pop()
    for parent, boundaries in inheritance_dict.items():
        for parent_bound, chimera_bound in boundaries.items():
            if not all(bound == 0 for bound in parent_bound):
                chi=sum(plddt_dict[chimera_label][slice(*chimera_bound)])
                native=sum(plddt_dict[parent][slice(*parent_bound)])
                raw_stability.append((chi - native) / native * 100)
    return np.mean(raw_stability)