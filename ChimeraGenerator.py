#!python
#
# Code 2023 by Jamel Simpson
"""Routines to generate chimeric sequences."""

from pathlib import Path
from json import load,dump


class chimeracls():
    pass


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
    # TODO make each sequence label malleable
    with open(file_name, 'w') as outfile:
        for (sequence, subunits, fasta_id) in list_of_sequence_subunits_label_tuples:
            for replicates in range(subunits):
                outfile.write(f'>{fasta_id}\n{sequence}\n')


def update_json(default_json, dilapidated_json):
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
    with open(str(Path(dilapidated_json).parent) + '/new' + str(Path(dilapidated_json).name), 'w') as f:
        dump(dilapidated, f, indent=4)


def change_value(dictionary,noted_key,new_value):
    for key,value in dictionary.items():
        if isinstance(value,dict):
            change_value(value,noted_key,new_value)
        elif key==noted_key:
            dictionary[noted_key]=new_value


def print_keys(dictionary):
    for key,value in dictionary.items():
        if isinstance(value,dict):
            print_keys(value)
        else:
            print(key)