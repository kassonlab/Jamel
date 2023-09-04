#!python
#
# Code 2023 by Jamel Simpson
"""Routines to generate chimeric sequences."""

from pathlib import Path
from json import load, dump
from os import path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

class chimeracls:
    pass

def general_attr_set(class_obj, dict_of_attrs):
    """All non-nested keys from the json are set as attributes to be called later"""
    for key, value in dict_of_attrs.items():
        setattr(class_obj, key, value)
    return class_obj


def sequence_splice(fasta_file, boundary_tuple, python_index='Yes'):
    """Takes a fasta sequence and returns the section of the sequence between indexes specified by the boundary one and two,
    as well as the sequence with the specified section replaced with a '-'.
    ABCDEFGH, boundary_one=0, boundary_two=3 Returns ABC and -DEFGH
    The residue at boundary_two is the first residue not included in the splice. If you're using alignment_finder, this is already accounted for."""
    if python_index == 'No':
        boundary_two = boundary_tuple[1] - 1
        boundary_one = boundary_tuple[0] - 1
    else:
        boundary_one = boundary_tuple[0]
        boundary_two = boundary_tuple[1]
    with open(fasta_file, "r") as fasta:
        sequence = ''.join(x for x in fasta if x[0] != '>' if x != '').strip().replace('\n', '')
    # The spliced region between the 2 specified boundaries is the first sequence
    # in the list followed by the sequence with the spliced region replace by a '-'
    spliced_region = ''.join(sequence[boundary_one:boundary_two])
    non_spliced_region = sequence.replace(spliced_region, '-')
    return spliced_region, non_spliced_region


def chimera_sequence_creation(section_being_spliced_in, marked_sequence, mark_in_the_sequence='-'):
    """Returns a completed sequence after replacing the '-' in a marked sequence with another sequence fragment"""
    chimera_sequence = marked_sequence.replace(mark_in_the_sequence, section_being_spliced_in)
    return chimera_sequence


def fasta_creation(file_name, list_of_sequence_subunits_label_tuples):
    """Creates a fasta file with the given file_name, and replicates the sequence within it the specified number of times
    to create a homo multimer if subunits is greater than 1."""
    sequences=[]
    with open(file_name, 'w') as outfile:

        for (sequence, subunits, fasta_id) in list_of_sequence_subunits_label_tuples:
            for replicates in range(subunits):
                sequences.append(SeqRecord(Seq(sequence), id=fasta_id,description=""))
        SeqIO.write(sequences, outfile, "fasta")


# TODO create a separate functionality file???
def update_json(default_json, dilapidated_json,overwrite=False):
    try:
        with open(default_json, 'r') as f:
            default = load(f)
        with open(dilapidated_json, 'r') as f:
            dilapidated = load(f)
    except FileNotFoundError:
        pass

    def merge_dict(default_dict, dilapidated_dict):
        for key, value in default_dict.items():
            if isinstance(value, dict) and key in dilapidated_dict and isinstance(dilapidated_dict[key], dict):
                merge_dict(value, dilapidated_dict[key])
            if key not in dilapidated_dict:
                dilapidated_dict[key] = value

    # TODO question marks ruin this

    def reduce_dict(default_dict, dilapidated_dict):
        for key, value in dilapidated_dict.copy().items():
            if key not in default_dict:
                del dilapidated_dict[key]
            if isinstance(value, dict) and key in default_dict and isinstance(default_dict[key], dict):
                reduce_dict(default_dict[key], value)

    merge_dict(default, dilapidated)
    reduce_dict(default, dilapidated)
    if overwrite:
        with open(dilapidated_json, 'w') as f:
            dump(dilapidated, f, indent=4)
    else:
        with open(str(Path(dilapidated_json).parent) + '/new' + str(Path(dilapidated_json).name), 'w') as f:
            dump(dilapidated, f, indent=4)


def change_json_value(json, noted_key='', new_value='',overwrite=False,find_n_replace_tuple=''):
    try:
        with open(json, 'r') as f:
            json_dict = load(f)
    except FileNotFoundError:
        pass

    def change(inner_dict):
        for key, value in inner_dict.items():
            if isinstance(value, dict):
                change(value)
            elif key == noted_key:
                inner_dict[key] = new_value
                break
    def find_n_replace(inner_dict):
        for key, value in inner_dict.items():
            if isinstance(value, dict):
                find_n_replace(value)
            elif isinstance(value, str) and find_n_replace_tuple[0] in value:
                print(value.replace(find_n_replace_tuple[0],find_n_replace_tuple[1]))
                inner_dict[key] = value.replace(find_n_replace_tuple[0],find_n_replace_tuple[1])
                break
    if find_n_replace_tuple:
        find_n_replace(json_dict)
    else:
        change(json_dict)
    if overwrite:
        with open(json, 'w') as f:
            dump(json_dict, f, indent=4)
    else:
        with open(str(Path(json).parent) + '/new' + str(Path(json).name), 'w') as f:
            dump(json_dict, f, indent=4)

def print_keys(dictionary,key_of_interest=''):
    dict_of_keys = {}

    def iterate_over_keys(inner_dictionary):
        for key, value in inner_dictionary.items():
            if isinstance(value, dict):
                iterate_over_keys(value)
            else:
                dict_of_keys[key]=value
    iterate_over_keys(dictionary)
    for key,value in sorted(list(dict_of_keys.items())):
        if key_of_interest and key_of_interest==key:
            print(key,':',value)
            break
    else:
        for key, value in sorted(list(dict_of_keys.items())):
            print(key)



def assign_file_attrs_to_chimeras(container):
    for chimera in container.chimeras:
        placeholder = container.naming_args.placeholder
        subunits = container.fasta_args.number_of_subunits
        alphafold_dir = container.naming_args.alphafold_outputs_dir
        chimera.monomer_stem = container.naming_args.monomer_naming_convention.replace(placeholder,
                                                                                       chimera.file_stem)
        chimera.chimera_stem = container.naming_args.chimera_naming_convention.replace(placeholder,
                                                                                       chimera.file_stem)
        chimera.chi_pdb = path.join(f'{alphafold_dir}{chimera.chimera_stem}', 'ranked_0.pdb')
        chimera.monomer_fasta = container.naming_args.fasta_directory + chimera.monomer_stem + container.naming_args.fasta_extension
        chimera.chimera_fasta = container.naming_args.fasta_directory + chimera.chimera_stem + container.naming_args.fasta_extension
        chimera.multimer_stem = container.naming_args.multimer_naming_convention.replace(placeholder,
                                                                                         chimera.file_stem)
        chimera.multimer_fasta = container.naming_args.fasta_directory + chimera.multimer_stem + container.naming_args.fasta_extension

        if subunits == 1:
            chimera.multimer_stem = container.naming_args.monomer_naming_convention.replace(placeholder,
                                                                                            chimera.file_stem)
            chimera.multimer_fasta = container.naming_args.fasta_directory + chimera.monomer_stem + container.naming_args.fasta_extension

        chimera.native_pdb = path.join(f'{alphafold_dir}{chimera.multimer_stem}', 'ranked_0.pdb')
        chimera.chi_pdb = path.join(f'{alphafold_dir}{chimera.chimera_stem}', 'ranked_0.pdb')
