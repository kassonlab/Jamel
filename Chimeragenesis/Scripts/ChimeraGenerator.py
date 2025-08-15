#!python
#
# Code 2023 by Jamel Simpson
"""Routines to generate chimeric sequences."""
import numpy as np
import pandas as pd

import AccessiontoAlignment
from pathlib import Path
from json import load, dump


def general_attr_set(class_obj, dict_of_attrs):
    """All non-nested keys from the json are set as attributes to be called later"""
    for key, value in dict_of_attrs.items():
        setattr(class_obj, key, value)
    return class_obj


def sequence_splice(sequence: str, splice_region: str, splice_marker: str = '#'):
    """Takes a protein sequence ABCDEFGH and splice_region ABC and Returns -DEFGH"""
    non_spliced_region = sequence.replace(splice_region, splice_marker)
    return non_spliced_region


def chimera_sequence_creation(section_being_spliced_in, marked_sequence, splice_marker='#'):
    """Returns a completed sequence after replacing the '-' in a marked sequence with another sequence fragment"""
    chimera_sequence = marked_sequence.replace(splice_marker, section_being_spliced_in)
    return chimera_sequence


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
        stem = Path(json_file).stem
        with open(json_file.replace(stem, f'new_{stem}'), 'w') as f:
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


def create_chimera_combinations(two_parent_aln_file: str, scanner_length, scanner_start=0, scanner_rate=1,
                                new_fasta_file=''):
    """Takes two sequence dictionary and creates all possible chimeras with given scanner parameters."""
    from ESM import SequenceDataframe
    aln_df = SequenceDataframe(two_parent_aln_file)
    parent1, parent2 = aln_df.index
    aln_df.add_value(parent1,'description',{parent1:None})
    aln_df.add_value(parent2, 'description', {parent2:None})

    def scanning_chimera_generator(base_label, partner_label):
        base_aln=aln_df.get_aln(base_label)
        splice_boundaries: list[tuple] = [(x, x + scanner_length) for x in
                                          range(scanner_start, len(base_aln)-scanner_length, scanner_rate)]
        base_splices=set(base_aln[slice(*boundary)].replace("-","") for boundary in splice_boundaries)
        for splice in base_splices:
            chimera_aln,inheritance_dict,parent_splice=AccessiontoAlignment.alignment_finder(splice, partner_label, base_label, two_parent_aln_file)
            parent_splice={key:str(value) for key,value in parent_splice.items()}
            aln_df.add_protein('_'.join([base_label,parent_splice["base_start"],parent_splice['base_end'],
                                         partner_label,parent_splice['partner_start'],parent_splice['partner_end']]),chimera_aln,inheritance_dict)

    scanning_chimera_generator(parent1, parent2)
    scanning_chimera_generator(parent2, parent1)
    if new_fasta_file:
        aln_df.dataframe_to_aln(new_fasta_file)
        aln_df.dataframe_to_multi_fa(Path(new_fasta_file).with_suffix('.mfa'))
    return aln_df


def create_combinations_no_aln(two_parent_fa_file: str, percentage_cutoff:tuple=(0.35,0.65),new_fasta_file=''):
    import sys
    sys.path.append('/scratch/jws6pq/Notebook/ESM/build/')
    import chi_cpp
    from ESM import SequenceDataframe
    from itertools import product
    aln_df = SequenceDataframe(two_parent_fa_file)
    parent1, parent2 = aln_df.index
    aln_df.add_value(parent1,'description',{parent1:None})
    aln_df.add_value(parent2, 'description', {parent2:None})

    # def single_cut_chimera_generator(base_label, partner_label):
    #     base_seq=aln_df.get_sequence(base_label)
    #     base_len=len(base_seq)
    #     base_cuts=range(*(round(base_len*fraction) for fraction in percentage_cutoff))
    #
    #     partner_seq = aln_df.get_sequence(partner_label)
    #     partner_len = len(partner_seq)
    #     partner_cuts = range(*(round(partner_len*fraction) for fraction in percentage_cutoff))
    #
    #     boundary_products=product(base_cuts,partner_cuts)
    #     for base_cut,partner_cut in boundary_products:
    #         base_splice=base_seq[0:base_cut]
    #         partner_splice=partner_seq[partner_cut:]
    #         chi_seq=base_splice+partner_splice
    #
    #         inheritance={base_label:{(0,len(base_splice)):(0,len(base_splice))},partner_label:{(partner_cut,None):(len(base_splice),None)}}
    #         aln_df.add_protein(f'{base_label}_0_{base_cut}_{partner_label}_{partner_cut}_{partner_len}',chi_seq,inheritance)
    # chi_cpp.single_cut_chimera_generator returns (chimera_labels:tuple,chimera_sequences:string,chimera_inheritance:dict)

    partner_label,partner_seq,partner_inh=chi_cpp.single_cut_chimera_generator(parent2,aln_df.get_sequence(parent2),parent1,aln_df.get_sequence(parent1), *percentage_cutoff)
    partner_info={'label':partner_label,"sequence":partner_seq,'aln_sequence':partner_seq,'description':partner_inh}
    base_label, base_seq, base_inh =chi_cpp.single_cut_chimera_generator(parent1,aln_df.get_sequence(parent1), parent2,aln_df.get_sequence(parent2),*percentage_cutoff)
    base_info={'label':base_label,"sequence":base_seq,'aln_sequence':base_seq,'description':base_inh}
    base_df=pd.DataFrame(base_info).set_index('label')
    partner_df=pd.DataFrame(partner_info).set_index('label')
    aln_df=aln_df._append(partner_df)
    aln_df=aln_df._append(base_df)
    aln_df.__class__=SequenceDataframe
    if new_fasta_file:
        aln_df.dataframe_to_aln(new_fasta_file)
        aln_df.dataframe_to_multi_fa(Path(new_fasta_file).with_suffix('.mfa'))
    return aln_df
