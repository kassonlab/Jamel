from ChimeraGenerator import fasta_creation
import Analysis
from AccessiontoAlignment import alignment_finder
from numpy import empty, savetxt
from os import system,path
from json import load
from sys import argv
from pathlib import Path
argument_json = argv[1]
with open(argument_json, 'rb') as jfile:
    argument_dict = load(jfile)["arguments"]
    operation_toggles = argument_dict['operation_toggles']
    fasta_arguments = argument_dict['fasta_arguments']
    analysis_arguments = argument_dict['analysis_arguments']
    alphafold_submission_args=argument_dict['alphafold_submission_args']
    check_arguments=argument_dict['check_arguments']
with open(fasta_arguments['alignment_file_name'], "r") as alignment:
    alignment = alignment.read().split('>')
    sequence_dictionary = {sequence.split('\n')[0]: ''.join(sequence.split('\n')[1:]) for sequence in alignment if
                           len(sequence) != 0}
scanner_length = fasta_arguments['scanner_length']
if scanner_length < 15:
    raise RuntimeError("scanner_length must be at 15 units")

naming_arguments = argument_dict['naming_arguments']
naming_convention = naming_arguments['naming_convention']
fasta_toggles = fasta_arguments['fasta_toggles']
ref_identifier = fasta_arguments['reference_identifier']
par_identifier = fasta_arguments['partner_identifier']
scanner_start = fasta_arguments['scanner_start']
scanner_end = scanner_start + scanner_length
scanner_movement_size = fasta_arguments['scanner_movement_size']
reference_sequence = sequence_dictionary[ref_identifier]
num_of_movements = fasta_arguments['num_of_movements']
no_gap_reference_sequence = ''.join(x for ind, x in enumerate(reference_sequence) if x.isalpha())
reference_length = len(no_gap_reference_sequence)
number_of_subunits = fasta_arguments['number_of_subunits'],
data_dict={}
data_dict['filename_stems']=data_dict['reference_boundaries']=data_dict['averaged_file_stems']=data_dict['reference_cuts']=data_dict['partner_cuts']=data_dict['partner_boundaries'] = ()
# TODO make a settings file taht saves settings for print
if fasta_toggles['Manually control number of scanner movements?'] == '':
    while scanner_end < reference_length:
        file_stem = naming_convention.replace(naming_arguments['reference_placeholder'],
                                              ref_identifier).replace(
            naming_arguments['partner_placeholder'], par_identifier).replace(
            naming_arguments['boundary_placeholder'], f'{scanner_start + 1}to{scanner_end}')
        data_dict['filename_stems']+=(file_stem,)
        averaged_stem = naming_arguments['averaged_plddt_naming_convention'].replace(
            naming_arguments['reference_placeholder'],
            ref_identifier).replace(
            naming_arguments['partner_placeholder'], par_identifier).replace(
            naming_arguments['boundary_placeholder'], f'{scanner_start + 1}to{scanner_end}')
        spliced_out = no_gap_reference_sequence[scanner_start:scanner_end]
        boundary_info = alignment_finder(fasta_arguments['alignment_file_name'], spliced_out, par_identifier,
                                         ref_identifier)
        data_dict['averaged_file_stems']+=(averaged_stem,)
        data_dict['reference_boundaries']+=(boundary_info[2],)
        data_dict['reference_cuts']+=(spliced_out,)
        data_dict['partner_cuts']+=(boundary_info[0],)
        data_dict['partner_boundaries']+=(boundary_info[1],)
        scanner_start += scanner_movement_size
        scanner_end += scanner_movement_size
else:
    count = 0
    while count <= num_of_movements:
        file_stem = naming_convention.replace(naming_arguments['reference_placeholder'],
                                              ref_identifier).replace(
            naming_arguments['partner_placeholder'], par_identifier).replace(
            naming_arguments['boundary_placeholder'], f'{scanner_start + 1}to{scanner_end}')
        data_dict['filename_stems']+=(file_stem,)
        averaged_stem = naming_arguments['averaged_plddt_naming_convention'].replace(
            naming_arguments['reference_placeholder'],
            ref_identifier).replace(
            naming_arguments['partner_placeholder'], par_identifier).replace(
            naming_arguments['boundary_placeholder'], f'{scanner_start + 1}to{scanner_end}')
        spliced_out = no_gap_reference_sequence[scanner_start:scanner_end]
        boundary_info = alignment_finder(fasta_arguments['alignment_file_name'], spliced_out, par_identifier,
                                         ref_identifier)
        data_dict['averaged_file_stems'] += (averaged_stem,)
        data_dict['reference_boundaries'] += (boundary_info[2],)
        data_dict['reference_cuts'] += (spliced_out,)
        data_dict['partner_cuts'] += (boundary_info[0],)
        data_dict['partner_boundaries'] += (boundary_info[1],)
        scanner_start += scanner_movement_size
        scanner_end += scanner_movement_size
        count += 1
fasta_list = ()
for stem in data_dict['filename_stems']:
    fasta_name = f'{naming_arguments["fasta_directory"]}{stem}{naming_arguments["fasta_file_extension"]}'
    fasta_list+=(fasta_name,)
# Matching python indexing for the indexing from the alignment
# with some amount of '-' and indexing in the regular sequence
if operation_toggles['run_fasta_operation?'] == '#':
    # Creating a regular sequence without '-'
    # Boundaries are given in python index
    for sequence_out, fasta,sequence_in in zip(data_dict['reference_cuts'], fasta_list,data_dict['partner_cuts']):
        spliced_sequence = no_gap_reference_sequence.replace(sequence_out, sequence_in)
        fasta_creation(fasta, spliced_sequence, number_of_subunits[0])
    if fasta_toggles['Make a list of created fasta files?'] == '#':
        fasta_list+=(fasta_arguments['reference_submission'],fasta_arguments['partner_submission'])
        with open(fasta_arguments['fasta_file_list_name'], 'w') as list_file:
            for fasta in fasta_list:
                list_file.write(f'{fasta}\n')
if operation_toggles['alphafold_submission'] == '#' or operation_toggles['check_submission'] == '#':
    submission_toggles = alphafold_submission_args['submission_toggles']
    check_toggles = check_arguments['check_toggles']
    fasta_list = list(fasta_list)
    output_directory = naming_arguments['alphafold_outputs_directory']
    if operation_toggles['check_submission'] == '#' or submission_toggles['use all fasta'] != '#':
        for index, fasta in enumerate(fasta_list):
            if path.exists(output_directory + Path(fasta).stem + '/ranking_debug.json'):
                del fasta_list[index]
    if check_toggles['create file of stragglers']=='#':
        with open(argument_dict['list_of_stragglers'], 'w') as run_list:
            for fasta in fasta_list:
                run_list.write(f'{fasta}')
    if fasta_list:
        proteins_per_slurm = alphafold_submission_args['proteins_per_slurm']
        template_slurm = alphafold_submission_args['template_slurm']
        alphafold_shell_script = alphafold_submission_args['alphafold_shell_script']
        naming_convention = naming_arguments['slurm_naming']
        placeholder = naming_arguments['slurm_placeholder']
        for slurm_index, file_index in enumerate(range(0, len(fasta_list), len(fasta_list) // proteins_per_slurm)):
            current_slurm = naming_convention.replace(placeholder, str(slurm_index))
            if submission_toggles['create_slurms'] == '#':
                system(f'cp {template_slurm} {current_slurm}')
                with open(current_slurm, 'a') as slurm_file:
                    slurm_file.write(
                        f'\n#SBATCH -o {alphafold_submission_args["slurm_output"].replace(placeholder, str(slurm_index))}\n'
                        f'#SBATCH -e {alphafold_submission_args["slurm_error"].replace(placeholder, str(slurm_index))}\n#Run program\n')
                    proteins_to_run = ','.join(fasta_list[file_index:file_index + proteins_per_slurm])
                    slurm_file.write(f'{alphafold_shell_script} {proteins_to_run} {output_directory}')
            if submission_toggles['sbatch slurms'] == '#':
                system(f'sbatch {current_slurm}')
# Analysis Portion
if operation_toggles['run_analysis_operation?'] == '#':
    analysis_data_dict={}
    analysis_toggles = analysis_arguments['analysis_toggles']
    analysis_data_dict['plddt_files'] = tuple(f'{naming_arguments["plddt_directory"]}{stem}{"plddt_file_extension"}' for stem in data_dict['filename_stems'])
    num_of_plddt=len(analysis_data_dict['plddt_files'])
    if analysis_toggles['make_plddts?'] == '#':
        analysis_data_dict['alphafold_folders'] = tuple(f"{naming_arguments['alphafold_outputs_directory']}{stem}/" for stem in
                             data_dict['filename_stems'])
        list(map(Analysis.generate_alphafold_files, analysis_data_dict['alphafold_folders'], analysis_data_dict['plddt_files']))
        reference_plddt = tuple(f"{naming_arguments['plddt_directory']}" \
                           f"{ref_identifier}{naming_arguments['plddt_file_extension']}")
        partner_plddt = tuple(f"{naming_arguments['plddt_directory']}{par_identifier}" \
                         f"{naming_arguments['plddt_file_extension']}")
        Analysis.generate_alphafold_files(f"{argument_dict['alphafold_outputs_directory']}"
                                          f"{ref_identifier}", reference_plddt[0])
        Analysis.generate_alphafold_files(f"{argument_dict['alphafold_outputs_directory']}{par_identifier}",
                                          partner_plddt[0])

        if number_of_subunits[0] > 1:
            analysis_data_dict['averaged_plddt'] = tuple(f'{naming_arguments["plddt_directory"]}{stem}{"plddt_file_extension"}' for stem in
                              data_dict['averaged_file_stems'])
            tuple(map(Analysis.averaging_multimer_plddt, analysis_data_dict['plddt_files'], analysis_data_dict['averaged_plddt'], number_of_subunits*num_of_plddt))
            averaged_reference = tuple(f"{naming_arguments['plddt_directory']}" \
                                  f"{naming_arguments['averaged_reference_plddt'].replace(naming_arguments['reference_placeholder'], ref_identifier)}{naming_arguments['plddt_file_extension']}")
            averaged_partner = tuple(f"{naming_arguments['plddt_directory']}" \
                                f"{naming_arguments['averaged_partner_plddt'].replace(naming_arguments['partner_placeholder'], ref_identifier)}{naming_arguments['plddt_file_extension']}")
            Analysis.averaging_multimer_plddt(f"{naming_arguments['plddt_directory']}"
                                              f"{ref_identifier}{naming_arguments['plddt_file_extension']}",
                                              averaged_reference[0],
                                              number_of_subunits[0])
            Analysis.averaging_multimer_plddt(f"{naming_arguments['plddt_directory']}"
                                              f"{par_identifier}{naming_arguments['plddt_file_extension']}",
                                              averaged_partner[0],
                                              number_of_subunits[0])

    if number_of_subunits[0] > 1:
        relative_stability = tuple(
            map(Analysis.average_relative_stability_full_chimera(averaged_partner, data_dict['partner_boundaries'],
                                                                 analysis_data_dict['averaged_plddt'], averaged_reference, data_dict['reference_cuts'],
                                                                 fasta_arguments['alignment_file_name'],
                                                                 ref_identifier)))
    else:
        relative_stability = tuple(
            map(Analysis.average_relative_stability_full_chimera(partner_plddt*num_of_plddt, data_dict['partner_boundaries'], analysis_data_dict['plddt_files'],
                                                                 reference_plddt*num_of_plddt, data_dict['reference_cuts'],
                                                                 fasta_arguments['alignment_file_name'],
                                                                 ref_identifier)))
    overall_chimera_stability = tuple(Analysis.overall_confidence(analysis_data_dict['plddt_files'])) * len(relative_stability)
    column_names = {key: value[1] for (key, value) in analysis_arguments['column_names'].items() if value[0] == '#'}
    data_array = empty((len(relative_stability) + 1, len(column_names)), dtype=object)
    for column_count, (corresponding_data, column_name) in enumerate(column_names.items()):
        exec(f'data_array[0, column_count], data_array[1:, column_count] = column_name, corresponding_data')
    savetxt(argument_dict['analysis_output_csv'], data_array, fmt=','.join('%s' for x in column_names), delimiter=",")

# lambda  chimera_class_tuple: map(func,chimera.attr1,chimera.attr2)
# from itertools import starmap
# map(func, kjmkj)
# starmap(func,((chi.attr1,chi.attr2)for chi in chi_list))
# using ** or * to explode the chimera tuple
from operator import attrgetter
from itertools import starmap
from ChimeraGenerator import fasta_creation, chimeracls
import Analysis
from AccessiontoAlignment import alignment_finder
from numpy import empty, savetxt
from os import system, path
from json import load
from sys import argv
from pathlib import Path

argument_json = argv[1]
with open(argument_json, 'rb') as jfile:
    argument_dict = load(jfile)["arguments"]
    operation_toggles = argument_dict['operation_toggles']
    fasta_arguments = argument_dict['fasta_arguments']
    analysis_arguments = argument_dict['analysis_arguments']
    alphafold_submission_args = argument_dict['alphafold_submission_args']
    check_arguments = argument_dict['check_arguments']
with open(fasta_arguments['alignment_file_name'], "r") as alignment:
    alignment = alignment.read().split('>')
    sequence_dictionary = {sequence.split('\n')[0]: ''.join(sequence.split('\n')[1:]) for sequence in alignment if
                           len(sequence) != 0}
scanner_length = fasta_arguments['scanner_length']
if scanner_length < 15:
    raise RuntimeError("scanner_length must be at 15 units")

naming_arguments = argument_dict['naming_arguments']
naming_convention = naming_arguments['naming_convention']
fasta_toggles = fasta_arguments['fasta_toggles']
ref_identifier = fasta_arguments['reference_identifier']
par_identifier = fasta_arguments['partner_identifier']
scanner_start = fasta_arguments['scanner_start']
scanner_end = scanner_start + scanner_length
scanner_movement_size = fasta_arguments['scanner_movement_size']
reference_sequence = sequence_dictionary[ref_identifier]
num_of_movements = fasta_arguments['num_of_movements']
no_gap_reference_sequence = ''.join(x for ind, x in enumerate(reference_sequence) if x.isalpha())
reference_length = len(no_gap_reference_sequence)
number_of_subunits = fasta_arguments['number_of_subunits']
chimera_container = ()
# TODO you need to use aceession alignment
# TODO tuples instead of list where possbile
# TODO make a settings file taht saves settings for print
if fasta_toggles['Manually control number of scanner movements?'] == '#':
    for index, end in enumerate(
            range(scanner_end, scanner_end + scanner_movement_size * num_of_movements, scanner_movement_size)):
        scanner_start = end - scanner_length
        exec(f'chimera_{index}=chimeracls()')
        chimera_container += globals()[f'chimera_{index}'],
        globals()[f'chimera_{index}'].filestem = naming_convention.replace(naming_arguments['reference_placeholder'],
                                                                           ref_identifier).replace(
            naming_arguments['partner_placeholder'], par_identifier).replace(
            naming_arguments['boundary_placeholder'], f'{scanner_start + 1}to{end}')
        globals()[f'chimera_{index}'].averaged_stem = naming_arguments['averaged_plddt_naming_convention'].replace(
            naming_arguments['reference_placeholder'],
            ref_identifier).replace(
            naming_arguments['partner_placeholder'], par_identifier).replace(
            naming_arguments['boundary_placeholder'], f'{scanner_start + 1}to{end}')
        spliced_out = no_gap_reference_sequence[scanner_start:end]
        boundary_info = alignment_finder(fasta_arguments['alignment_file_name'], spliced_out, par_identifier,
                                         ref_identifier)
        globals()[f'chimera_{index}'].reference_boundaries = boundary_info[2]
        globals()[f'chimera_{index}'].reference_cuts = spliced_out
        globals()[f'chimera_{index}'].partner_cuts = boundary_info[0]
        globals()[f'chimera_{index}'].partner_boundaries = boundary_info[1]
else:
    for index, end in enumerate(range(scanner_end, reference_length, scanner_movement_size)):
        scanner_start = end - scanner_length
        exec(f'chimera_{index}=chimeracls()')
        chimera_container += globals()[f'chimera_{index}'],
        globals()[f'chimera_{index}'].file_stem = naming_convention.replace(naming_arguments['reference_placeholder'],
                                                                            ref_identifier).replace(
            naming_arguments['partner_placeholder'], par_identifier).replace(
            naming_arguments['boundary_placeholder'], f'{scanner_start + 1}to{end}')
        globals()[f'chimera_{index}'].averaged_stem = naming_arguments['averaged_plddt_naming_convention'].replace(
            naming_arguments['reference_placeholder'],
            ref_identifier).replace(
            naming_arguments['partner_placeholder'], par_identifier).replace(
            naming_arguments['boundary_placeholder'], f'{scanner_start + 1}to{end}')
        spliced_out = no_gap_reference_sequence[scanner_start:end]
        boundary_info = alignment_finder(fasta_arguments['alignment_file_name'], spliced_out, par_identifier,
                                         ref_identifier)
        globals()[f'chimera_{index}'].reference_boundaries = boundary_info[2]
        globals()[f'chimera_{index}'].reference_cuts = spliced_out
        globals()[f'chimera_{index}'].partner_cuts = boundary_info[0]
        globals()[f'chimera_{index}'].partner_boundaries = boundary_info[1]
num_of_chis = len(chimera_container)
for chimera in chimera_container:
    chimera.fasta_name = f'{naming_arguments["fasta_directory"]}{chimera.filestem}{naming_arguments["fasta_file_extension"]}'
if operation_toggles['run_fasta_operation?'] == '#':
    for chimera in chimera_container:
        spliced_sequence = no_gap_reference_sequence.replace(chimera.reference_cuts, chimera.partner_cuts)
        fasta_creation(chimera.fasta_name, spliced_sequence, number_of_subunits)
    if fasta_toggles['Make a list of created fasta files?'] == '#':
        with open(fasta_arguments['fasta_file_list_name'], 'w') as list_file:
            for chimera in chimera_container:
                list_file.write(f'{chimera.fasta_name}\n')
            list_file.write(fasta_arguments['reference_submission'])
            list_file.write(fasta_arguments['partner_submission'])
if operation_toggles['run_analysis_operation?'] == '#':
    analysis_toggles = analysis_arguments['analysis_toggles']
    reference_plddt = f"{naming_arguments['plddt_directory']}" \
                      f"{ref_identifier}{naming_arguments['plddt_file_extension']}"
    partner_plddt = f"{naming_arguments['plddt_directory']}{par_identifier}" \
                    f"{naming_arguments['plddt_file_extension']}"
    averaged_reference = f"{naming_arguments['plddt_directory']}" \
                         f"{naming_arguments['averaged_reference_plddt'].replace(naming_arguments['reference_placeholder'], ref_identifier)}{naming_arguments['plddt_file_extension']}"
    averaged_partner = f"{naming_arguments['plddt_directory']}" \
                       f"{naming_arguments['averaged_partner_plddt'].replace(naming_arguments['partner_placeholder'], ref_identifier)}{naming_arguments['plddt_file_extension']}"
    for chimera in chimera_container:
        chimera.plddt = f'{naming_arguments["plddt_directory"]}{chimera.file_stem}{"plddt_file_extension"}'
        if number_of_subunits > 1:
            chimera.avgplddt = f'{naming_arguments["plddt_directory"]}{chimera.averaged_stem}{"plddt_file_extension"}'

    if analysis_toggles['make_plddts?'] == '#':
        Analysis.generate_alphafold_files(f"{argument_dict['alphafold_outputs_directory']}"
                                          f"{ref_identifier}", reference_plddt)
        Analysis.generate_alphafold_files(f"{argument_dict['alphafold_outputs_directory']}{par_identifier}",
                                          partner_plddt)
        if number_of_subunits > 1:
            Analysis.averaging_multimer_plddt(f"{naming_arguments['plddt_directory']}"
                                              f"{ref_identifier}{naming_arguments['plddt_file_extension']}",
                                              averaged_reference,
                                              number_of_subunits)
            Analysis.averaging_multimer_plddt(f"{naming_arguments['plddt_directory']}"
                                              f"{par_identifier}{naming_arguments['plddt_file_extension']}",
                                              averaged_partner,
                                              number_of_subunits)
            for chimera in chimera_container:
                Analysis.generate_alphafold_files(
                    f"{naming_arguments['alphafold_outputs_directory']}{chimera.file_stem}/", chimera.plddt)
                if number_of_subunits > 1:
                    Analysis.averaging_multimer_plddt(chimera.plddt, chimera.avgplddt, number_of_subunits)
    for chimera in chimera_container:
        if number_of_subunits[0] > 1:
            chimera.rel_stability = Analysis.average_relative_stability_full_chimera(averaged_partner,
                                                                                     chimera.partner_boundaries,
                                                                                     chimera.avgplddt,
                                                                                     averaged_reference,
                                                                                     chimera.reference_cuts,
                                                                                     fasta_arguments[
                                                                                         'alignment_file_name'],
                                                                                     ref_identifier)
        else:
            chimera.rel_stability = Analysis.average_relative_stability_full_chimera(partner_plddt,
                                                                                     chimera.partner_boundaries,
                                                                                     chimera.plddt,
                                                                                     reference_plddt,
                                                                                     chimera.reference_cuts,
                                                                                     fasta_arguments[
                                                                                         'alignment_file_name'],
                                                                                     ref_identifier)
    data_dict = {
        'overall_chimera_stability': tuple(Analysis.overall_confidence(chimera.plddt) for chimera in chimera_container),
        'relative_stability': tuple(chimera.plddt for chimera in chimera_container),
        'filename_stems': tuple(chimera.file_stem for chimera in chimera_container)}
    column_names = {key: value[1] for (key, value) in analysis_arguments['column_names'].items() if value[0] == '#'}
    data_array = empty((num_of_chis + 1, len(column_names)), dtype=object)
    for column_count, (corresponding_data, column_name) in enumerate(column_names.items()):
        data_array[0, column_count], data_array[1:, column_count] = column_name, data_dict[corresponding_data]
    savetxt(argument_dict['analysis_output_csv'], data_array, fmt=','.join('%s' for x in column_names), delimiter=",")
    del data_dict