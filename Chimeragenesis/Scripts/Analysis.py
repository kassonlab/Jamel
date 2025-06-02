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
import pandas as pd
from matplotlib import pyplot as plt


def scatterplot_for_df(df:pd.DataFrame, x, y, attr_mods:dict=None, scatter_args:dict=None):
    if attr_mods is None:
        attr_mods = {plt.title: f'{x} vs. {y}'}
    if scatter_args is None:
        scatter_args = {}
    df.plot.scatter(x=x, y=y,**scatter_args)
    if attr_mods:
        for attr,mod in attr_mods.items():
            attr(mod)
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

def embedding_heatmap(embedding_df:pd.DataFrame,attr_mods:dict=None, scatter_args:dict=None,new_csv=''):
    """Creates a scatter plot where every chimera in the embedding df is plotted by the length that protein A
    donated to its n-terminus in the x-axis and the length that protein B donated to its c-terminus in the y-axis"""
    if scatter_args is None:
        scatter_args={'c': 'dot_product', 'colormap': 'viridis', 's': 50, 'xlabel': 'N-term Length',
         'ylabel': 'C-term Length'}
    if attr_mods is None:
        attr_mods = {plt.title: f'Dot Product Heat Map'}
    embedding_df=embedding_df.assign(base_len=lambda series: [[int(end) - int(start) for start, end in re.findall(r'_(\d+)_(\d+)', index)][0] for index in series.index],
                                     partner_len=lambda series: [[int(end) - int(start) for start, end in re.findall(r'_(\d+)_(\d+)', index)][1] for index in series.index],
                                     chimera_len= lambda series: series['sequence'].apply(lambda element: len(element)))
    scatterplot_for_df(embedding_df, 'base_len', 'partner_len',scatter_args=scatter_args,attr_mods=attr_mods)
    if new_csv:
        embedding_df.to_csv(new_csv,index_label='label')

def embedding_dp_vs_length(embedding_df:pd.DataFrame,attr_mods:dict=None, scatter_args:dict=None):
    if attr_mods is None:
        attr_mods = {plt.title: f'Chimera Length vs. Dot Product'}
    embedding_df = embedding_df.assign(
        chimera_len=lambda series: series['sequence'].apply(lambda element: len(element)))
    scatterplot_for_df(embedding_df, 'chimera_len', 'dot_product',scatter_args=scatter_args,attr_mods=attr_mods)
def select_for_af_from_embedding_df(embedding_df:pd.DataFrame,amount_from_top:int,total_labels_wanted:int,file_for_selected_labels:str,parents:list):
    """Based on cutoffs puts the corresponding chimera labels in a text file for alphafold submission"""
    embedding_df = embedding_df.sort_values(by='dot_product', ascending=False)
    search_interval=embedding_df.__len__()//(total_labels_wanted-amount_from_top)
    with open(file_for_selected_labels, "w") as file:
        print(*(embedding_df['dot_product'].iloc[:amount_from_top].index.to_list() + embedding_df['dot_product'].iloc[amount_from_top::search_interval].index.to_list()+parents),
              sep='\n',file=file)
def process_raw_embedding_csv(embedding_csv:str,parents:list):
    embedding_df=pd.read_csv(embedding_csv,index_col='label')
    embedding_df = embedding_df[~embedding_df.index.duplicated(keep='first')]
    embedding_df.to_csv(embedding_csv,index_label='label')
    for parent in parents:
        parent_df = embedding_df[embedding_df.index.str.startswith(parent)]
        parent_df = parent_df.drop(index=parent)
        parent_df = parent_df.sort_values(by='dot_product')
        embedding_heatmap(parent_df,scatter_args={'c':'dot_product','colormap':'viridis','s':25,'xlabel':f'{parent} N-term Length','ylabel':f'{parent} C-term Length'})
        embedding_dp_vs_length(parent_df,scatter_args={'s':25,'xlabel':'Chimera Length','ylabel':'Dot Product'},attr_mods={plt.title:f'{parent} Chis Len v DP'})

def graph_plddt_per_residue(pdb_list:list[str],attr_mods:dict=None, scatter_args:dict=None):
    if attr_mods is None:
        attr_mods = {plt.title: f'Residue # vs. Confidence',plt.xlabel: 'Residue #',plt.ylabel:"Confidence Score"}
    if scatter_args is None:
        scatter_args = {}
    for pdb in pdb_list:
        for seq,plddt in get_plddt_dict_from_pdb(pdb).items():
            plt.plot(range(len(seq)),plddt,label=Path(pdb).stem,**scatter_args)

    if attr_mods:
        for attr,mod in attr_mods.items():
            attr(mod)

    plt.legend()
    plt.show()
