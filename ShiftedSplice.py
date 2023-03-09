from sys import argv
from json import load
from ChimeraGenerator import fasta_creation
import Analysis
from AccessiontoAlignment import alignment_finder
from numpy import empty, savetxt

argument_json = argv[1]
with open(argument_json, 'rb') as jfile:
    argument_dict = load(jfile)["arguments"]
    operation_toggles = argument_dict['operation_toggles']
    fasta_arguments = argument_dict['fasta_arguments']
    analysis_arguments = argument_dict['analysis_arguments']
with open(fasta_arguments['alignment_file_name'], "r") as alignment:
    alignment = alignment.read().split('>')
    # Splitting up the sequences names and sequences into a dictionary
    sequence_dictionary = {sequence.split('\n')[0]: ''.join(sequence.split('\n')[1:]) for sequence in alignment if
                           len(sequence) != 0}
scanner_length = fasta_arguments['scanner_length']
if scanner_length < 15:
    raise RuntimeError("scanner_length must be at 15 units")
naming_settings = argument_dict['naming_arguments']
naming_convention = naming_settings['naming_convention']
fasta_toggles = fasta_arguments['fasta_toggles']
# define number of movements
ref_identifier = fasta_arguments['reference_identifier']
par_identifier = fasta_arguments['partner_identifier']
scanner_start = fasta_arguments['scanner_start']
scanner_end = scanner_start + scanner_length
scanner_movement_size = fasta_arguments['scanner_movement_size']
reference_sequence = sequence_dictionary[ref_identifier]
num_of_movements = fasta_arguments['num_of_movements'][1]
no_gap_reference_sequence = ''.join(x for ind, x in enumerate(reference_sequence) if x.isalpha())
reference_length = len(no_gap_reference_sequence)
number_of_subunits = fasta_arguments['number_of_subunits']
filename_stems = []
reference_boundaries = []
averaged_file_stems = []
reference_cuts = []
partner_cuts = []
partner_boundaries = []  # you need to use aceession alignment
if fasta_toggles['Manually control number of scanner movements?'] == '':
    while scanner_end < reference_length:
        file_stem = naming_convention.replace(fasta_arguments['reference_placeholder'],
                                              ref_identifier).replace(
            fasta_arguments['partner_placeholder'], par_identifier).replace(
            fasta_arguments['boundary_placeholder'], f'{scanner_start + 1}to{scanner_end}')
        filename_stems.append(file_stem)
        averaged_stem = naming_settings['averaged_plddt_naming_convention'].replace(
            fasta_arguments['reference_placeholder'],
            ref_identifier).replace(
            fasta_arguments['partner_placeholder'], par_identifier).replace(
            fasta_arguments['boundary_placeholder'], f'{scanner_start + 1}to{scanner_end}')
        spliced_out = no_gap_reference_sequence[scanner_start:scanner_end]
        boundary_info = alignment_finder(fasta_arguments['alignment_file_name'], spliced_out, par_identifier,
                                         ref_identifier)
        averaged_file_stems.append(averaged_stem)
        reference_boundaries.append(boundary_info[2])
        reference_cuts.append(spliced_out)
        partner_cuts.append(boundary_info[0])
        partner_boundaries.append(boundary_info[1])
        scanner_start += scanner_movement_size
        scanner_end += scanner_movement_size
else:
    count = 0
    while count <= num_of_movements:
        file_stem = naming_convention.replace(naming_settings['reference_placeholder'],
                                              ref_identifier).replace(
            naming_settings['partner_placeholder'], par_identifier).replace(
            naming_settings['boundary_placeholder'], f'{scanner_start + 1}to{scanner_end}')
        filename_stems.append(file_stem)
        averaged_stem = naming_settings['averaged_plddt_naming_convention'].replace(
            fasta_arguments['reference_placeholder'],
            ref_identifier).replace(
            fasta_arguments['partner_placeholder'], par_identifier).replace(
            fasta_arguments['boundary_placeholder'], f'{scanner_start + 1}to{scanner_end}')
        spliced_out = no_gap_reference_sequence[scanner_start:scanner_end]
        boundary_info = alignment_finder(fasta_arguments['alignment_file_name'], spliced_out, par_identifier,
                                         ref_identifier)
        averaged_file_stems.append(averaged_stem)
        reference_boundaries.append(boundary_info[2])
        reference_cuts.append(spliced_out)
        partner_cuts.append(boundary_info[0])
        partner_boundaries.append(boundary_info[1])
        scanner_start += scanner_movement_size
        scanner_end += scanner_movement_size
        count += 1
# Matching python indexing for the indexing from the alignment
# with some amount of '-' and indexing in the regular sequence
if operation_toggles['run_fasta_operation?'] == '#':
    # Creating a regular sequence without '-'
    # Boundaries are given in python index
    fasta_list = []
    for sequence_out, stem, sequence_in in zip(reference_cuts, filename_stems, partner_cuts):
        spliced_sequence = no_gap_reference_sequence.replace(sequence_out, sequence_in)
        fasta_name = f'{naming_settings["fasta_directory"]}{stem}{naming_settings["fasta_file_extension"]}'
        fasta_list.append(fasta_name)
        fasta_creation(fasta_name, spliced_sequence, number_of_subunits)
    if fasta_toggles['Make a list of created fasta files?'] == '#':
        fasta_list.append()
        with open(fasta_arguments['fasta_file_list_name'], 'w') as list_file:
            for fasta in fasta_list:
                list_file.write(f'{fasta}\n')

# Analysis Portion
if operation_toggles['run_analysis_operation?'] == '#':
    analysis_toggles = analysis_arguments['analysis_toggles']
    plddt_files = [f'{naming_settings["plddt_directory"]}{stem}{"plddt_file_extension"}' for stem in filename_stems]
    if analysis_toggles['make_plddts?'] == '#':
        alphafold_folders = [f"{argument_dict['alphafold_outputs_directory']}{stem}/" for stem in
                             filename_stems]
        list(map(Analysis.generate_alphafold_files, alphafold_folders, plddt_files))
        reference_plddt = [f"{naming_settings['plddt_directory']}" \
                           f"{ref_identifier}{naming_settings['plddt_file_extension']}" for x in filename_stems]
        partner_plddt = [f"{naming_settings['plddt_directory']}{par_identifier}" \
                         f"{naming_settings['plddt_file_extension']}" for x in filename_stems]
        Analysis.generate_alphafold_files(f"{argument_dict['alphafold_outputs_directory']}"
                                          f"{ref_identifier}", reference_plddt[0])
        Analysis.generate_alphafold_files(f"{argument_dict['alphafold_outputs_directory']}{par_identifier}",
                                          partner_plddt[0])

        if number_of_subunits[0] > 1:
            averaged_plddt = [f'{naming_settings["plddt_directory"]}{stem}{"plddt_file_extension"}' for stem in
                              averaged_file_stems]
            list(map(Analysis.averaging_multimer_plddt, plddt_files, averaged_plddt, number_of_subunits))
            averaged_reference = [f"{naming_settings['plddt_directory']}" \
                                  f"{naming_settings['averaged_reference_plddt'].replace(fasta_arguments['reference_placeholder'], ref_identifier)}{naming_settings['plddt_file_extension']}"
                                  for x in filename_stems]
            averaged_partner = [f"{naming_settings['plddt_directory']}" \
                                f"{naming_settings['averaged_partner_plddt'].replace(fasta_arguments['partner_placeholder'], ref_identifier)}{naming_settings['plddt_file_extension']}"
                                for x in filename_stems]
            Analysis.averaging_multimer_plddt(f"{naming_settings['plddt_directory']}"
                                              f"{ref_identifier}{naming_settings['plddt_file_extension']}",
                                              averaged_reference[0],
                                              number_of_subunits[0])
            Analysis.averaging_multimer_plddt(f"{naming_settings['plddt_directory']}"
                                              f"{par_identifier}{naming_settings['plddt_file_extension']}",
                                              averaged_partner[0],
                                              number_of_subunits[0])

    if number_of_subunits[0] > 1:
        relative_stability = list(
            map(Analysis.average_relative_stability_full_chimera(averaged_partner, partner_boundaries,
                                                                 averaged_plddt, averaged_reference, reference_cuts,
                                                                 fasta_arguments['alignment_file_name'],
                                                                 ref_identifier)))
    else:
        relative_stability = list(
            map(Analysis.average_relative_stability_full_chimera(partner_plddt, partner_boundaries, plddt_files,
                                                                 reference_plddt, reference_cuts,
                                                                 fasta_arguments['alignment_file_name'],
                                                                 ref_identifier)))
    overall_chimera_stability = [Analysis.overall_confidence(plddt_files)] * len(relative_stability)
    column_names = {key: value[1] for (key, value) in analysis_arguments['column_names'].items() if value[0] == '#'}
    exec(f'')
    data_array = empty((len(relative_stability) + 1, len(column_names)), dtype=object)
    for column_count, (corresponding_data, column_name) in enumerate(column_names.items()):
        exec(f'data_array[0, column_count], data_array[1:, column_count] = column_name, corresponding_data')

    savetxt(argument_dict['analysis_output_csv'], data_array, fmt=','.join('%s' for x in column_names), delimiter=",")
