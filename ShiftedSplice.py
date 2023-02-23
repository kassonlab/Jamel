from sys import argv
from json import load
from ChimeraGenerator import fasta_creation

argument_json = argv[1]
with open(argument_json, 'rb') as jfile:
    argument_dict = load(jfile)["arguments"]
with open(argument_dict['alignment_file_name'], "r") as alignment:
    alignment = alignment.read().split('>')
    # Splitting up the sequences names and sequences into a dictionary
    sequence_dictionary = {sequence.split('\n')[0]: ''.join(sequence.split('\n')[1:]) for sequence in alignment if
                           len(sequence) != 0}
num_of_movements = argument_dict['num_of_movements'][1]
scanner_length = argument_dict['scanner_length']
scanner_start = argument_dict['scanner_start']
scanner_movement_size = argument_dict['scanner_movement_size']
scanner_end = scanner_start + scanner_length
reference_sequence = sequence_dictionary[argument_dict['reference_identifier']]
comparison_sequence = list(sequence_dictionary[argument_dict['partner_identifier']])
number_of_subunits = argument_dict['number_of_subunits']
# Matching python indexing for the indexing from the alignment
# with some amount of '-' and indexing in the regular sequence
reference_sequence_indexing = [ind for ind, x in enumerate(reference_sequence) if x.isalpha()]
# Creating a regular sequence without '-'
no_gap_reference_sequence = [x for ind, x in enumerate(reference_sequence) if x.isalpha()]
# Boundaries are given in python index
fasta_list=[]
if argument_dict['num_of_movements'][0] != '':
    while scanner_end < len(reference_sequence_indexing):
        spliced_sequence = no_gap_reference_sequence.copy()
        alignment_index = (reference_sequence_indexing[scanner_start], reference_sequence_indexing[scanner_end])
        spliced_sequence[scanner_start:scanner_end] = comparison_sequence[alignment_index[0]:alignment_index[1]]
        spliced_sequence = ''.join(spliced_sequence).replace('-', '')
        fasta_name = argument_dict['naming_convention'].replace(argument_dict['reference_placeholder'],
                                                                argument_dict['reference_identifier']).replace(
            argument_dict['partner_placeholder'], argument_dict['partner_identifier']).replace(
            argument_dict['boundary_placeholder'], f'{scanner_start + 1}to{scanner_end}')
        fasta_list.append(fasta_name)
        fasta_creation(fasta_name, spliced_sequence, number_of_subunits)
        scanner_start += scanner_movement_size
        scanner_end += scanner_movement_size
else:
    count = 0
    while count <= num_of_movements:
        spliced_sequence = no_gap_reference_sequence.copy()
        alignment_index = (reference_sequence_indexing[scanner_start], reference_sequence_indexing[scanner_end])
        spliced_sequence[scanner_start:scanner_end] = comparison_sequence[alignment_index[0]:alignment_index[1]]
        spliced_sequence = ''.join(spliced_sequence).replace('-', '')
        fasta_name = argument_dict['naming_convention'].replace(argument_dict['reference_placeholder'],
                                                                argument_dict['reference_identifier']).replace(
            argument_dict['partner_placeholder'], argument_dict['partner_identifier']).replace(
            argument_dict['boundary_placeholder'], f'{scanner_start + 1}to{scanner_end}')
        fasta_list.append(fasta_name)
        fasta_creation(fasta_name, spliced_sequence, number_of_subunits)
        scanner_start += scanner_movement_size
        scanner_end += scanner_movement_size
        count += 1
with open(argument_dict['fasta_file_list_name'][1], 'w') as list_file:
    for fasta in fasta_list:
        list_file.write(f'{fasta}\n')
