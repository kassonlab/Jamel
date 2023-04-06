
from json import load
from os import system, path
from pathlib import Path
from sys import argv
from time import perf_counter
from numpy import empty, savetxt
import Analysis
from setup import create_alphafold_slurm
from AccessiontoAlignment import alignment_finder
from ChimeraGenerator import fasta_creation, chimeracls
from itertools import product


# This class creates a container for the chimeras and settings generated for each protein and associated argument json file

class ShiftedChimeraContainer:
    def __init__(self, shifted_json_file):
        with open(shifted_json_file, 'rb') as jfile:
            self.argument_dict = load(jfile)["arguments"]
        self.operation_toggles = self.argument_dict['operation_toggles']
        self.fasta_arguments = self.argument_dict['fasta_arguments']
        self.scanner_length = self.fasta_arguments['scanner_length']
        if self.scanner_length < 10:
            raise RuntimeError("scanner_length must be at 10 units")
        self.analysis_arguments = self.argument_dict['analysis_arguments']
        self.alphafold_submission_args = self.argument_dict['alphafold_submission_args']
        self.check_arguments = self.argument_dict['check_arguments']
        self.ref_identifier = self.fasta_arguments['reference_identifier']
        with open(self.fasta_arguments['alignment_file_name'], "r") as alignment:
            alignment = alignment.read().split('>')
            sequence_dictionary = {sequence.split('\n')[0]: ''.join(sequence.split('\n')[1:]) for sequence in alignment
                                   if
                                   len(sequence) != 0}
            self.reference_sequence = sequence_dictionary[self.ref_identifier]
        self.naming_arguments = self.argument_dict['naming_arguments']
        self.naming_convention = self.naming_arguments['naming_convention']
        self.fasta_toggles = self.fasta_arguments['fasta_toggles']
        self.par_identifier = self.fasta_arguments['partner_identifier']
        self.scanner_start = self.fasta_arguments['scanner_start']
        self.scanner_end = self.scanner_start + self.scanner_length
        self.scanner_movement_size = self.fasta_arguments['scanner_movement_size']
        self.num_of_movements = self.fasta_arguments['num_of_movements']
        self.no_gap_reference_sequence = ''.join(x for ind, x in enumerate(self.reference_sequence) if x.isalpha())
        self.reference_length = len(self.no_gap_reference_sequence)
        self.number_of_subunits = self.fasta_arguments['number_of_subunits']
        self.chimeras = ()

    def add_chimera(self, chimera):
        self.chimeras += (chimera,)


one = perf_counter()
arg_jsons = argv[1:]
container_of_containers = ()
for index, jsn in enumerate(arg_jsons):
    container_of_containers += (ShiftedChimeraContainer(jsn),)
# Heteromers are created by submitting more than one json argument, the first submitted json will be used for operation settings and toggles, not pertinent to creating the chimera sequence
prime_container = container_of_containers[0]
# TODO make a settings file taht saves settings for print
for container in container_of_containers:
    if container.fasta_toggles['Manually control number of scanner movements?'] != '#':
        for index, end in enumerate(
                range(container.scanner_end, container.reference_length, container.scanner_movement_size)):
            scanner_start = end - container.scanner_length
            container.add_chimera(chimeracls())

            # This if statement is checking for the boundary placeholder intended for naming processes.
            # If it's left blank, the boundaries from this group of chimeras/protein will be left out of naming schemes.
            # This is very useful for cases where you dont want one of the proteins in the heteromer to be a chimera and therefore has no splice boundaries
            if container.naming_arguments['boundary_placeholder']:
                container.chimeras[index].no_gap_boundaries = f'{scanner_start + 1}to{end}'

            spliced_out = container.no_gap_reference_sequence[scanner_start:end]
            boundary_info = alignment_finder(container.fasta_arguments['alignment_file_name'], spliced_out,
                                             container.par_identifier,
                                             container.ref_identifier)
            container.chimeras[index].reference_boundaries = boundary_info[2]
            container.chimeras[index].reference_cuts = spliced_out
            container.chimeras[index].partner_cuts = boundary_info[0]
            container.chimeras[index].partner_boundaries = boundary_info[1]
            container.chimeras[index].container_id = container
    else:
        for index, end in enumerate(
                range(container.scanner_end,
                      container.scanner_end + container.scanner_movement_size * container.num_of_movements,
                      container.scanner_movement_size)):
            scanner_start = end - container.scanner_length
            container.add_chimera(chimeracls())
            spliced_out = container.no_gap_reference_sequence[scanner_start:end]

            # This if statement is checking for the boundary placeholder intended for naming processes.
            # If it's left blank, the boundaries from this group of chimeras/protein will be left out of naming schemes.
            # This is very useful for cases where you dont want one of the proteins in the heteromer to be a chimera and therefore has no splice boundaries
            if container.naming_arguments['boundary_placeholder']:
                container.chimeras[index].no_gap_boundaries = f'{scanner_start + 1}to{end}'
            boundary_info = alignment_finder(container.fasta_arguments['alignment_file_name'], spliced_out,
                                             container.par_identifier,
                                             container.ref_identifier)
            container.chimeras[index].reference_boundaries = boundary_info[2]
            container.chimeras[index].reference_cuts = spliced_out
            container.chimeras[index].partner_cuts = boundary_info[0]
            container.chimeras[index].partner_boundaries = boundary_info[1]
            container.chimeras[index].container_id = container
    num_of_chis = len(container.chimeras)

# This if/else statement creates every combination between all protein list generated by each json file,
if prime_container.fasta_toggles["Make pair or combo heteromers"] == "combo" and len(container_of_containers)>1:
    groupings = tuple(product(*(container.chimeras for container in container_of_containers)))
# otherwise all list will be paired one to one to one etc. or single file list will be created for homomeric proteins
else:
    groupings = tuple(zip(*(container.chimeras for container in container_of_containers)))
for chimera_group in groupings:
    joined_boundary = '_'.join(
        chimera.no_gap_boundaries for chimera in chimera_group if hasattr(chimera, 'no_gap_boundaries'))
    chimera_group[0].file_stem = prime_container.naming_convention.replace(
        prime_container.naming_arguments['reference_placeholder'],
        prime_container.ref_identifier).replace(
        prime_container.naming_arguments['partner_placeholder'], prime_container.par_identifier).replace(
        prime_container.naming_arguments['boundary_placeholder'], joined_boundary)
    chimera_group[0].averaged_stem = prime_container.naming_arguments['averaged_plddt_naming_convention'].replace(
        prime_container.naming_arguments['reference_placeholder'],
        prime_container.ref_identifier).replace(
        prime_container.naming_arguments['partner_placeholder'], prime_container.par_identifier).replace(
        prime_container.naming_arguments['boundary_placeholder'], joined_boundary)
    chimera_group[0].fasta_name = path.join(prime_container.naming_arguments["fasta_directory"],
                                            f'{chimera_group[0].file_stem}{prime_container.naming_arguments["fasta_file_extension"]}')
    for chimera in chimera_group:
        chimera.chi_sequence=chimera.container_id.no_gap_reference_sequence.replace(chimera.reference_cuts,
                                                               chimera.partner_cuts)
# TODO make the fasta labels more accurate???
if prime_container.operation_toggles['run_fasta_operation?'] == '#':
    for chimera_group in groupings:
        spliced_sequences = []
        for chimera in chimera_group:
            spliced_sequences.append((chimera.chi_sequence,
                                      chimera.container_id.number_of_subunits))
        fasta_creation(chimera_group[0].fasta_name, spliced_sequences)
if prime_container.fasta_toggles['Make a list of created fasta files?'] == '#':
    with open(prime_container.fasta_arguments['fasta_file_list_name'], 'w') as list_file:
        list_file.write('\n'.join(chimera_group[0].fasta_name for chimera_group in groupings))
        list_file.write(
            '\n' + prime_container.fasta_arguments['reference_submission'] + '\n' + prime_container.fasta_arguments[
                'partner_submission'])
if prime_container.operation_toggles['alphafold_submission'] == '#' or prime_container.operation_toggles[
    'check_submission'] == '#':
    submission_toggles = prime_container.alphafold_submission_args['submission_toggles']
    check_toggles = prime_container.check_arguments['check_toggles']
    output_directory = prime_container.naming_arguments['alphafold_outputs_directory']
    fastas = [chimera_group[0].fasta_name for chimera_group in groupings] + [
        container.fasta_arguments['reference_submission'] for container in container_of_containers] + [
                 container.fasta_arguments['partner_submission'] for container in container_of_containers]
    fasta_to_run = ()
    # This operation will check the completion of an alphafold prediction and return a tuple with the fastas with incomplete predictions,
    # this tuple can either be turn into list in a file or ran by the next alphafold_submission operation
    if prime_container.operation_toggles['check_submission'] == '#':
        for fasta in fastas:
            if not path.exists(output_directory + Path(fasta).stem + '/ranking_debug.json'):
                fasta_to_run += (fasta,)
        if check_toggles['create file of stragglers'] == '#':
            with open(prime_container.check_arguments['list_of_stragglers'], 'w') as run_list:
                for fasta in fasta_to_run:
                    run_list.write('\n'.join(fasta for fasta in fasta_to_run))
        if not fasta_to_run:
            prime_container.operation_toggles['alphafold_submission'] = ''
    # If check_submission operation isn't run
    #TODO allow list of stragglers to be run, as an option along with check submission
    if prime_container.operation_toggles['alphafold_submission'] == '#':
        proteins_per_slurm = prime_container.alphafold_submission_args['proteins_per_slurm']
        template_slurm = prime_container.alphafold_submission_args['template_slurm']
        alphafold_shell_script = prime_container.alphafold_submission_args['alphafold_shell_script']
        naming_convention = prime_container.naming_arguments['slurm_naming']
        placeholder = prime_container.naming_arguments['slurm_placeholder']
        if not fasta_to_run:
            fasta_to_run = fastas
        for slurm_index, file_index in enumerate(range(0, len(fasta_to_run), proteins_per_slurm)):
            current_slurm = naming_convention.replace(placeholder, str(slurm_index))

            if submission_toggles['create_slurms'] == '#':
                create_alphafold_slurm(fasta_to_run[file_index:file_index + proteins_per_slurm], current_slurm,
                                       template_slurm,
                                       prime_container.alphafold_submission_args["slurm_output"].replace(placeholder,
                                                                                                         str(slurm_index)),
                                       prime_container.alphafold_submission_args["slurm_error"].replace(placeholder,
                                                                                                        str(slurm_index)),
                                       alphafold_shell_script, output_directory)
            if submission_toggles['sbatch slurms'] == '#':
                system(f'sbatch {current_slurm}')
if prime_container.operation_toggles['run_analysis_operation?'] == '#':
    analysis_toggles = prime_container.analysis_arguments['analysis_toggles']
    pdbs = [path.join(f'{prime_container.naming_arguments["alphafold_outputs_directory"]}{Path(container.fasta_arguments["reference_submission"]).stem}','ranked_0.pdb') for container in container_of_containers] + [
        path.join(f'{prime_container.naming_arguments["alphafold_outputs_directory"]}{Path(container.fasta_arguments["partner_submission"]).stem}','ranked_0.pdb') for container in container_of_containers]
    for chimera_group in groupings:
        pdb=path.join(f'{prime_container.naming_arguments["alphafold_outputs_directory"]}{chimera_group[0].file_stem}', 'ranked_0.pdb')
        pdbs.append(pdb)
        plddt_dict=Analysis.get_plddt_from_pdb(pdb)
        for chimera in chimera_group:
            if chimera.chi_sequence in plddt_dict:
                chimera.plddt=plddt_dict[chimera.chi_sequence]
            else:
                print(chimera.file_stem, pdb)
    # reference_plddt = f"{prime_container.naming_arguments['plddt_directory']}" \
    #                   f"{Path(fasta_arguments['reference_submission']).stem}{prime_container.naming_arguments['plddt_file_extension']}"
    # partner_plddt = f"{prime_container.naming_arguments['plddt_directory']}{Path(fasta_arguments['partner_submission']).stem}" \
    #                 f"{prime_container.naming_arguments['plddt_file_extension']}"
    # averaged_reference = f"{prime_container.naming_arguments['plddt_directory']}" \
    #                      f"{prime_container.naming_arguments['averaged_reference_plddt'].replace(naming_arguments['reference_placeholder'], Path(fasta_arguments['reference_submission']).stem)}{prime_container.naming_arguments['plddt_file_extension']}"
    # averaged_partner = f"{prime_container.naming_arguments['plddt_directory']}" \
    #                    f"{naming_arguments['averaged_partner_plddt'].replace(naming_arguments['partner_placeholder'], Path(fasta_arguments['partner_submission']).stem)}{prime_container.naming_arguments['plddt_file_extension']}"
    # for chimera in chimera_container:
    #     chimera.plddt = f'{prime_container.naming_arguments["plddt_directory"]}{chimera.file_stem}{prime_container.naming_arguments["plddt_file_extension"]}'
    #     if number_of_subunits > 1:
    #         chimera.avgplddt = f'{prime_container.naming_arguments["plddt_directory"]}{chimera.averaged_stem}{prime_container.naming_arguments["plddt_file_extension"]}'
    #
    # if analysis_toggles['make_plddts?'] == '#':
    #     Analysis.generate_alphafold_files(f"{prime_container.naming_arguments['alphafold_outputs_directory']}"
    #                                       f"{Path(fasta_arguments['reference_submission']).stem}/", reference_plddt)
    #     Analysis.generate_alphafold_files(f"{prime_container.naming_arguments['alphafold_outputs_directory']}{Path(fasta_arguments['partner_submission']).stem}/",
    #                                       partner_plddt)
    #     if number_of_subunits > 1:
    #         Analysis.averaging_multimer_plddt(reference_plddt,
    #                                           averaged_reference,
    #                                           number_of_subunits)
    #         Analysis.averaging_multimer_plddt(partner_plddt,
    #                                           averaged_partner,
    #                                           number_of_subunits)
    #         for chimera in chimera_container:
    #             Analysis.generate_alphafold_files(
    #                 f"{prime_container.naming_arguments['alphafold_outputs_directory']}{chimera.file_stem}/", chimera.plddt)
    #             if number_of_subunits > 1:
    #                 Analysis.averaging_multimer_plddt(chimera.plddt, chimera.avgplddt, number_of_subunits)

    # for chimera in chimera_container:
    #     if number_of_subunits > 1:
    #         chimera.rel_stability = Analysis.average_relative_stability_full_chimera(averaged_partner,
    #                                                                                  chimera.partner_boundaries,
    #                                                                                  chimera.avgplddt,
    #                                                                                  averaged_reference,
    #                                                                                  chimera.reference_cuts,
    #                                                                                  prime_container.fasta_arguments[
    #                                                                                      'alignment_file_name'],
    #                                                                                  ref_identifier)
    #     else:
    #         chimera.rel_stability = Analysis.average_relative_stability_full_chimera(partner_plddt,
    #                                                                                  chimera.partner_boundaries,
    #                                                                                  chimera.plddt,
    #                                                                                  reference_plddt,
    #                                                                                  chimera.reference_cuts,
    #                                                                                  prime_container.fasta_arguments[
    #                                                                                      'alignment_file_name'],
    #                                                                                  ref_identifier)
    # data_dict = {
    #     'overall_chimera_stability': tuple(Analysis.overall_confidence(chimera.plddt) for chimera in chimera_container),
    #     'relative_stability': tuple(chimera.plddt for chimera in chimera_container),
    #     'filename_stems': tuple(chimera.file_stem for chimera in chimera_container)}
    # column_names = {key: value[1] for (key, value) in analysis_arguments['column_names'].items() if value[0] == '#'}
    # data_array = empty((num_of_chis + 1, len(column_names)), dtype=object)
    # for column_count, (corresponding_data, column_name) in enumerate(column_names.items()):
    #     data_array[0, column_count], data_array[1:, column_count] = column_name, data_dict[corresponding_data]
    # savetxt(analysis_arguments['analysis_output_file'], data_array, fmt=','.join('%s' for x in column_names), delimiter=",")
    # del data_dict
two = perf_counter()
print(two - one)
