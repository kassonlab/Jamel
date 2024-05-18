#! python
#
# Code 2023 by Jamel Simpson


from pickle import load as p_load
from json import load as j_load
from os import system, path, strerror, listdir, makedirs
from shutil import copy
from numpy import savetxt, empty,zeros
from errno import ENOENT
from Bio import PDB
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

def convert_array_to_file(arr,delimiter,file_name):
    savetxt(file_name, arr, fmt=delimiter.join(('%s' for x in arr)))
def get_plddt_dict_from_pdb(pdb_file):
    pdb = PDB.PDBParser().get_structure('pdb', pdb_file)[0]
    homomeric = defaultdict(tuple)
    chains = tuple(chain for chain in pdb)
    for chain in chains:
        sequence = ''.join(PDB.Polypeptide.three_to_one(resi.get_resname()) for resi in chain)
        plddt = tuple(resi['CA'].bfactor for resi in chain)
        homomeric[sequence] += (plddt,)
    homomeric = {sequence: tuple(sum(scores) / len(scores) for scores in zip(*homomers)) for sequence, homomers in
                 homomeric.items()}
    return homomeric


def get_plddt_file_from_pdb(pdb_file, new_plddt_file):
    pdb = PDB.PDBParser().get_structure('pdb', pdb_file)[0]
    homomeric = defaultdict(tuple)
    chains = tuple(chain for chain in pdb)
    for chain in chains:
        sequence = ''.join(PDB.Polypeptide.three_to_one(resi.get_resname()) for resi in chain)
        plddt = tuple(resi['CA'].bfactor for resi in chain)
        homomeric[sequence] += (plddt,)
    homomeric = {sequence: tuple(str(round(sum(scores) / len(scores), 2)) for scores in zip(*homomers)) for
                 sequence, homomers in homomeric.items()}
    with open(new_plddt_file, 'w') as new_plddt:
        new_plddt.write(
            '\n'.join('>{0}\n{1}'.format(sequence, "\n".join(plddt)) for sequence, plddt in homomeric.items()))


def skeletonize_alphafold_folder(alphafold_dir, storage_dir):
    rank_file = path.join(alphafold_dir, 'ranking_debug.json')
    makedirs(storage_dir, exist_ok=True)
    new_files = []
    dir_files = [file for file in listdir(alphafold_dir) if file.startswith('ranked')]
    try:
        copy(rank_file, path.join(storage_dir, 'ranking_debug.json'))
        for pdb_file in dir_files:
            new_pdb = path.join(storage_dir, pdb_file)
            new_plddt = path.join(storage_dir, str(Path(pdb_file).stem) + '.plddt')
            new_files.append(new_pdb)
            new_files.append(new_plddt)
            copy(path.join(alphafold_dir, pdb_file), new_pdb)
            get_plddt_file_from_pdb(path.join(alphafold_dir, pdb_file), new_plddt)
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


def generate_alphafold_files(alphafold_folder, new_plddt='', new_pdb=''):
    """Creates a text file containing the plddt values of the highest_rank_model extracted from alphafold's result pkl file
    and renames the ranked_0.pdb file and places it in the desired directory."""
    # Checking to see if ranking_debug.json exists. This file is the last to be output by alphafold and is a check that
    # the pkl file you want to extract from exists, as well as to avoid errors
    ranking_file = path.join(alphafold_folder, 'ranking_debug.json')
    try:
        if new_pdb:
            # The highest ranked structure is copied with a new name and directory
            copy(path.join(alphafold_folder, 'ranked_0.pdb'), new_pdb)
        if new_plddt:
            # ranking_debug is also useful for determining which result pkl file is the highest ranked. The model result pkl files are
            # numbered by the order they are created and not their overall confidence score. The information about their rank by
            # confidence score is found in ranking_debug.json
            with open(ranking_file, 'r') as jfile:
                highest_rank_model = j_load(jfile)['order'][0]
            with open(path.join(alphafold_folder, f'result_{highest_rank_model}.pkl'), 'rb') as pfile:
                # The plddt scores are put into a column in a text file named by new_plddt
                savetxt(new_plddt, p_load(pfile)['plddt'], fmt='%s', delimiter=' ')
    except FileNotFoundError:
        raise FileNotFoundError(ENOENT, strerror(ENOENT), ranking_file)


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


def overall_confidence(plddt_tuple):
    """Returns the average confidence score from a protein's plddt file."""
    average_plddt = sum(plddt_tuple) / len(plddt_tuple)
    return average_plddt

def turn_plddt_dict_into_tuples(plddt_dict):
    plddt_list_of_tuples = []
    for seq, plddt in plddt_dict.items():
        plddt_list_of_tuples.append((seq, plddt))
    return plddt_list_of_tuples


def get_chimera_boundaries(chimera_seq, seq_spliced_into_ref):
    if chimera_seq.find(seq_spliced_into_ref) != -1:
        splice_start = chimera_seq.find(seq_spliced_into_ref)
    else:
        print('broken')
        return
    splice_end = splice_start + len(seq_spliced_into_ref)
    # Then those boundaries are compared against the very beginning and end of the proteins, by introducing them into a set
    # to get of redundancy if the sequence_of_interest boundaries contain the beginning or end

    boundaries = tuple({0, splice_start, splice_end, len(chimera_seq)})
    boundaries = sorted(boundaries)
    chimera_boundaries = []
    for index in range(len(boundaries) - 1):
        if (boundaries[index], boundaries[index + 1]) == (splice_start, splice_end):
            chimera_boundaries.append(('native', boundaries[index], boundaries[index + 1]))
        else:
            chimera_boundaries.append(('ref', boundaries[index], boundaries[index + 1]))
    return chimera_boundaries


def relative_stability(native_plddt, native_boundary_tuple, chimera_plddt, chimera_boundary_tuple):
    """Returns the relative percent difference between the two equally sized sections of plddt scores that are outlined with
    native_boundary_tuple and chimera_boundary_tuple. relative_difference=(compared value-reference value)/reference value * 100
    Native scores are assumed to be the reference value in this formula for relative difference"""
    # Pulling the plddt values as floats that start at native_boundary_tuple[0] and chimera_boundary_tuple[0], and end at
    # native_boundary_tuple[1] and chimera_boundary_tuple[1] but dont include index [1] scores.
    native_score = native_plddt[native_boundary_tuple[0]:native_boundary_tuple[1]]
    chimera_score = chimera_plddt[chimera_boundary_tuple[0]:chimera_boundary_tuple[1]]
    # Recording the length of the residue scores for averaging purposes later
    splice_length = len(chimera_score)
    relative_difference = sum((chimera - native) / native * 100 for native, chimera in zip(native_score, chimera_score))
    return relative_difference, splice_length


# def turn_rank_matrix_into_dict(rank_matrix_file):


# TODO compare sequences before lloking at relative stability
def revamped_rs(native_plddt_dict, chimera_plddt_dict, reference_plddt_dict, seq_spliced_into_ref):
    raw_stability = 0
    native_seq = turn_plddt_dict_into_tuples(native_plddt_dict)[0][0]
    native_plddt = turn_plddt_dict_into_tuples(native_plddt_dict)[0][1]
    chi_seq = turn_plddt_dict_into_tuples(chimera_plddt_dict)[0][0]
    chi_plddt = turn_plddt_dict_into_tuples(chimera_plddt_dict)[0][1]
    reference_seq = turn_plddt_dict_into_tuples(reference_plddt_dict)[0][0]
    reference_plddt = turn_plddt_dict_into_tuples(reference_plddt_dict)[0][1]
    chimera_boundaries = get_chimera_boundaries(chi_seq, seq_spliced_into_ref)
    for boundary in chimera_boundaries:
        seq_chunk = chi_seq[boundary[1]:boundary[2]]
        if boundary[0] == 'native':
            start = native_seq.find(seq_chunk)
            end = start + len(seq_chunk)
            native_chunk = native_seq[start:end]
            if native_chunk == seq_chunk:
                raw_stability += relative_stability(native_plddt, (start, end), chi_plddt, boundary[1:])[0]
            else:
                print('sequences arent equal')
                break
        elif boundary[0] == 'ref':
            start = reference_seq.find(seq_chunk)
            end = start + len(seq_chunk)
            ref_chunk = reference_seq[start:end]
            if ref_chunk == seq_chunk:
                raw_stability += relative_stability(reference_plddt, (start, end), chi_plddt, boundary[1:])[0]
            else:
                print('sequences arent equal')
                break
        else:
            print('sequences arent equal')
            break
    return raw_stability / len(chi_seq)

def averaging_multimer_plddt(plddt_file, new_plddt_file, subunits):
    """This function takes a plddt and averages the scores
    for each residue position across the number of subunints specified"""
    # Using list comprehension to turn the plddt file into a list of floats
    with open(plddt_file, 'r') as infile:
        multimer_plddt = tuple(float(score) for score in infile.readlines())
    # Calculating the length a subunits to have for step size when iterating through the list later
    monomer_length = int(len(multimer_plddt) / int(subunits))
    # using list comprehension to step through each the residue position of each subunit and
    # collect their scores, average them and return them to the new list
    averaged_scores = tuple(sum(multimer_plddt[residue_index::monomer_length]) / subunits
                            for residue_index in range(monomer_length))
    # creating a file to input the averaged scores
    with open(new_plddt_file, 'w') as new_plddt:
        new_plddt.write('\n'.join(str(score) for score in averaged_scores))






