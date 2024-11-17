#!python
#
# Code 2023 by Jamel Simpson
"""Routines to generate chimeric sequences."""
import AccessiontoAlignment
from pathlib import Path
from json import load, dump
from os import path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO





def general_attr_set(class_obj, dict_of_attrs):
    """All non-nested keys from the json are set as attributes to be called later"""
    for key, value in dict_of_attrs.items():
        setattr(class_obj, key, value)
    return class_obj


def sequence_splice(sequence:str, splice_region:str, splice_marker: str = '#'):
    """Takes a protein sequence and Returns ABC and -DEFGH
    The residue at boundary_two is the first residue not included in the splice. If you're using alignment_finder, this is already accounted for."""

    # The spliced region between the 2 specified boundaries is the first sequence
    # in the list followed by the sequence with the spliced region replace by a str marker like '#'
    non_spliced_region = sequence.replace(splice_region, splice_marker)
    return non_spliced_region


def chimera_sequence_creation(section_being_spliced_in, marked_sequence, splice_marker='#'):
    """Returns a completed sequence after replacing the '-' in a marked sequence with another sequence fragment"""
    chimera_sequence = marked_sequence.replace(splice_marker, section_being_spliced_in)
    return chimera_sequence

def get_chimera(base_aln_seq,seq_to_splice:str):
    marked_base=sequence_splice(AccessiontoAlignment.no_gap_sequence_from_alignment(base_aln_seq),seq_to_splice)
    splice_region=AccessiontoAlignment.alignment_finder()



def update_json(default_json, dilapidated_json, overwrite=False):
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


def change_json_value(json_file, noted_key='', new_value='', overwrite=False, find_n_replace_tuple=''):
    """Iterates through json keys and replaces their value accordingly"""
    try:
        with open(json_file, 'r') as f:
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
                print(value.replace(find_n_replace_tuple[0], find_n_replace_tuple[1]))
                inner_dict[key] = value.replace(find_n_replace_tuple[0], find_n_replace_tuple[1])
                break

    if find_n_replace_tuple:
        find_n_replace(json_dict)
    else:
        change(json_dict)
    if overwrite:
        with open(json_file, 'w') as f:
            dump(json_dict, f, indent=4)
    else:
        stem=Path(json_file).stem
        with open(json_file.replace(stem,f'new_{stem}'), 'w') as f:
            dump(json_dict, f, indent=4)


def print_keys(dictionary, key_of_interest=''):
    dict_of_keys = {}

    def iterate_over_keys(inner_dictionary):
        for key, value in inner_dictionary.items():
            if isinstance(value, dict):
                iterate_over_keys(value)
            else:
                dict_of_keys[key] = value

    iterate_over_keys(dictionary)
    for key, value in sorted(list(dict_of_keys.items())):
        if key_of_interest and key_of_interest == key:
            print(key, ':', value)
            break
    else:
        for key, value in sorted(list(dict_of_keys.items())):
            print(key)


def create_chimera_combinations(two_seq_dict: dict, scanner_length, scanner_start=0, scanner_rate=1, new_fasta_file='') -> dict[str, str]:
    """Takes two sequence dictionary and creates all possible chimeras with given splice length.
    Returns dict[parent1#parent2#splice1#splice2, chimera_sequence]"""
    from AccessiontoAlignment import dictionary_to_fasta
    new_chimera_dict = {label.replace('-',''):seq for label,seq in two_seq_dict.items()}

    def scanning_chimera_generator(base_sequence, partner_seq, base_label, partner_label):
        splice_boundaries: list[tuple] = [(x, x + scanner_length) for x in
                                          range(scanner_start, len(base_sequence) - scanner_length, scanner_rate) if
                                          x + scanner_length < len(base_sequence)]
        for boundary in splice_boundaries:
            # This is the sequence to be spliced into by the partner seq, the residues between a given boundary have
            # been cut and replaced with a marker character '-' for the aligned sequence from the partner to replace it
            marked_base_sequence = sequence_splice(base_sequence, boundary)[1]
            partner_replacement = sequence_splice(partner_seq, boundary)[0]
            chimeric_seq = chimera_sequence_creation(partner_replacement, marked_base_sequence)
            new_chimera_dict['-'.join((base_label,partner_label,str(boundary[0]),str(boundary[1])))] = chimeric_seq
        return new_chimera_dict

    seq1, seq2 = new_chimera_dict.values()
    seq1_label, seq2_label = new_chimera_dict.keys()
    scanning_chimera_generator(seq1, seq2, seq1_label, seq2_label)
    scanning_chimera_generator(seq2, seq1, seq2_label, seq1_label)
    if new_fasta_file:
        dictionary_to_fasta(new_chimera_dict, new_fasta_file)
    return new_chimera_dict
