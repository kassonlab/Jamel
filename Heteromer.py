# need to be able to store boundaries and sequences  for every unique seqeucen
# unique alignments for every unqie sequence
# Most likely previous shifted script eith a for loop for every unique sequence
# from itertools import groupby
# with open('E:\Jamel\\ranked_0.pdb', 'r') as pdb:
#     pdb=[tuple(line.split()) for line in pdb.readlines() if line.split()[0] not in ('TER','END')]
#     for index,line in enumerate(pdb):
#         print(index,line[10])
#     print(list(groupby(pdb,lambda x:x[5])))
# TODO create different proteins by creating a new object everytime you find a one, create a sequence out of the amino acids match the sequences given and the sequences found, use plddt from pdb for calculations
# TODO create a json file for every unqiue protein
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
prime_container = container_of_containers[0]
# TODO make a settings file taht saves settings for print
for container in container_of_containers:
    if container.fasta_toggles['Manually control number of scanner movements?'] != '#':
        for index, end in enumerate(
                range(container.scanner_end, container.reference_length, container.scanner_movement_size)):
            scanner_start = end - container.scanner_length
            container.add_chimera(chimeracls())
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
if prime_container.fasta_toggles["Make pair or combo heteromers"] == "combo":
    groupings = tuple(product(*(container.chimeras for container in container_of_containers)))
else:
    groupings = zip(*(container.chimeras for container in container_of_containers))
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
# TODO make the fasta labels more accurate???
if prime_container.operation_toggles['run_fasta_operation?'] == '#':
    for chimera_group in groupings:
        spliced_sequences = []
        for chimera in chimera_group:
            spliced_sequences.append((chimera.container_id.no_gap_reference_sequence.replace(chimera.reference_cuts,
                                                                                             chimera.partner_cuts),
                                      chimera.container_id.number_of_subunits))
        fasta_creation(chimera_group[0].fasta_name, spliced_sequences)
if prime_container.fasta_toggles['Make a list of created fasta files?'] == '#':
    with open(prime_container.fasta_arguments['fasta_file_list_name'], 'w') as list_file:
        list_file.write('\n'.join(chimera_group[0].fasta_name for chimera_group in groupings))
        list_file.write(
            '\n' + prime_container.fasta_arguments['reference_submission'] + '\n' + prime_container.fasta_arguments[
                'partner_submission'])
if prime_container.operation_toggles['alphafold_submission'] == '#' or prime_container.operation_toggles['check_submission'] == '#':
    submission_toggles = prime_container.alphafold_submission_args['submission_toggles']
    check_toggles = prime_container.check_arguments['check_toggles']
    output_directory = prime_container.naming_arguments['alphafold_outputs_directory']
    fastas = [chimera_group[0].fasta_name for chimera_group in groupings] + [
        prime_container.fasta_arguments['reference_submission']] + [
                 prime_container.fasta_arguments['partner_submission']]
    fasta_to_run = ()
    if prime_container.operation_toggles['check_submission'] == '#':
        for fasta in fastas:
            if not path.exists(output_directory + Path(fasta).stem + '/ranking_debug.json'):
                print(Path(fasta).stem)
                fasta_to_run += (fasta,)
        if check_toggles['create file of stragglers'] == '#':
            with open(prime_container.check_arguments['list_of_stragglers'], 'w') as run_list:
                for fasta in fasta_to_run:
                    run_list.write('\n'.join(fasta for fasta in fasta_to_run))
        if not fasta_to_run:
            prime_container.operation_toggles['alphafold_submission'] = ''
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
# if operation_toggles['run_analysis_operation?'] == '#':
#     analysis_toggles = analysis_arguments['analysis_toggles']
#     reference_plddt = f"{naming_arguments['plddt_directory']}" \
#                       f"{Path(fasta_arguments['reference_submission']).stem}{naming_arguments['plddt_file_extension']}"
#     partner_plddt = f"{naming_arguments['plddt_directory']}{Path(fasta_arguments['partner_submission']).stem}" \
#                     f"{naming_arguments['plddt_file_extension']}"
#     averaged_reference = f"{naming_arguments['plddt_directory']}" \
#                          f"{naming_arguments['averaged_reference_plddt'].replace(naming_arguments['reference_placeholder'], Path(fasta_arguments['reference_submission']).stem)}{naming_arguments['plddt_file_extension']}"
#     averaged_partner = f"{naming_arguments['plddt_directory']}" \
#                        f"{naming_arguments['averaged_partner_plddt'].replace(naming_arguments['partner_placeholder'], Path(fasta_arguments['partner_submission']).stem)}{naming_arguments['plddt_file_extension']}"
#     for chimera in chimera_container:
#         chimera.plddt = f'{naming_arguments["plddt_directory"]}{chimera.file_stem}{naming_arguments["plddt_file_extension"]}'
#         if number_of_subunits > 1:
#             chimera.avgplddt = f'{naming_arguments["plddt_directory"]}{chimera.averaged_stem}{naming_arguments["plddt_file_extension"]}'
#
#     if analysis_toggles['make_plddts?'] == '#':
#         Analysis.generate_alphafold_files(f"{naming_arguments['alphafold_outputs_directory']}"
#                                           f"{Path(fasta_arguments['reference_submission']).stem}/", reference_plddt)
#         Analysis.generate_alphafold_files(f"{naming_arguments['alphafold_outputs_directory']}{Path(fasta_arguments['partner_submission']).stem}/",
#                                           partner_plddt)
#         if number_of_subunits > 1:
#             Analysis.averaging_multimer_plddt(reference_plddt,
#                                               averaged_reference,
#                                               number_of_subunits)
#             Analysis.averaging_multimer_plddt(partner_plddt,
#                                               averaged_partner,
#                                               number_of_subunits)
#             for chimera in chimera_container:
#                 Analysis.generate_alphafold_files(
#                     f"{naming_arguments['alphafold_outputs_directory']}{chimera.file_stem}/", chimera.plddt)
#                 if number_of_subunits > 1:
#                     Analysis.averaging_multimer_plddt(chimera.plddt, chimera.avgplddt, number_of_subunits)
#     for chimera in chimera_container:
#         if number_of_subunits > 1:
#             chimera.rel_stability = Analysis.average_relative_stability_full_chimera(averaged_partner,
#                                                                                      chimera.partner_boundaries,
#                                                                                      chimera.avgplddt,
#                                                                                      averaged_reference,
#                                                                                      chimera.reference_cuts,
#                                                                                      fasta_arguments[
#                                                                                          'alignment_file_name'],
#                                                                                      ref_identifier)
#         else:
#             chimera.rel_stability = Analysis.average_relative_stability_full_chimera(partner_plddt,
#                                                                                      chimera.partner_boundaries,
#                                                                                      chimera.plddt,
#                                                                                      reference_plddt,
#                                                                                      chimera.reference_cuts,
#                                                                                      fasta_arguments[
#                                                                                          'alignment_file_name'],
#                                                                                      ref_identifier)
#     data_dict = {
#         'overall_chimera_stability': tuple(Analysis.overall_confidence(chimera.plddt) for chimera in chimera_container),
#         'relative_stability': tuple(chimera.plddt for chimera in chimera_container),
#         'filename_stems': tuple(chimera.file_stem for chimera in chimera_container)}
#     column_names = {key: value[1] for (key, value) in analysis_arguments['column_names'].items() if value[0] == '#'}
#     data_array = empty((num_of_chis + 1, len(column_names)), dtype=object)
#     for column_count, (corresponding_data, column_name) in enumerate(column_names.items()):
#         data_array[0, column_count], data_array[1:, column_count] = column_name, data_dict[corresponding_data]
#     savetxt(analysis_arguments['analysis_output_file'], data_array, fmt=','.join('%s' for x in column_names), delimiter=",")
#     del data_dict
two = perf_counter()
print(two - one)
