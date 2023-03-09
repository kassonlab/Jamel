import ChimeraGenerator
from sys import argv
from json import load
import AccessiontoAlignment
from numpy import savetxt, empty
import Analysis
from pathlib import Path

argument_json = argv[1]
protein_info = argv[2]
# ADD Filesort and alphafold submission and gromacs simulation

with open(argument_json, 'rb') as jfile:
    argument_dict = load(jfile)["arguments"]
    fasta_arguments = argument_dict['fasta_arguments']
    operation_toggles = argument_dict['operation_toggles']
    prefix_protein = argument_dict['Which protein is first, constant or variant?']
    naming_arguments = argument_dict['naming_arguments']
    analysis_arguments = argument_dict['analysis_arguments']
with open(protein_info, 'r') as info_list:
    info_list = info_list.readlines()
    accession_number = [x.split()[0] for x in info_list]
    protein_names = [x.split()[1] for x in info_list]
constant_sequence_of_interest = fasta_arguments['constant_sequence_of_interest']
variant_sequence_of_interest = fasta_arguments['variant_sequence_of_interest']
with open(constant_sequence_of_interest, 'r') as fasta:
    constant_sequence_of_interest = ''.join([x for x in fasta if x[0] != '>' if x != '']).strip().replace('\n', '')
with open(variant_sequence_of_interest, 'r') as fasta:
    variant_sequence_of_interest = ''.join([x for x in fasta if x[0] != '>' if x != '']).strip().replace('\n', '')

constant_fasta = fasta_arguments['constant_fasta']
variant_fasta = fasta_arguments['variant_fasta']
placeholder = naming_arguments['placeholder']
constant_fasta_identifier = fasta_arguments['constant_fasta_identifier']
variant_fasta_identifier = [fasta_arguments['variant_fasta_identifier'] for name in protein_names]
msa = [fasta_arguments['msa_file_name'] for name in protein_names]
constant_sequence_of_interest = [constant_sequence_of_interest for name in protein_names]
variant_sequence_of_interest = [variant_sequence_of_interest for name in protein_names]
number_of_subunits = [fasta_arguments['number_of_subunits'] for x in protein_names]

if operation_toggles['run_fasta_operation?'] == '#':
    email = [fasta_arguments['email_for_accession'] for name in protein_names]
    monomer_fastas = [naming_arguments['fasta_directory'] + naming_arguments['monomer_naming_convention'].replace(
        placeholder, name) + naming_arguments['fasta_file_extension'] for name in protein_names]
    chimera_fastas = [naming_arguments['fasta_directory'] + naming_arguments['chimera_naming_convention'].replace(
        placeholder, name) + naming_arguments['fasta_file_extension'] for name in protein_names]
    fasta_toggles = fasta_arguments['fasta_toggles']
    if number_of_subunits[0] == 1:
        list(map(AccessiontoAlignment.accession_to_fasta, monomer_fastas, accession_number, email, number_of_subunits))
    else:
        multimer_fastas = [naming_arguments['fasta_directory'] + naming_arguments['multimer_naming_convention'].replace(
            placeholder, name) + naming_arguments['fasta_file_extension'] for name in protein_names]
        list(map(AccessiontoAlignment.accession_to_fasta, monomer_fastas, accession_number, email, number_of_subunits,
                 multimer_fastas))

    if fasta_toggles['Create an alignment?'] == '#':
        AccessiontoAlignment.multiple_sequence_alignment(monomer_fastas, fasta_arguments['msa_fasta'], msa[0],
                                                         variant_fasta, fasta_arguments['muscle_command_for_msa'])
    variant_sequences = [x[0] for x in
                         map(AccessiontoAlignment.alignment_finder, msa, variant_sequence_of_interest, protein_names,
                             variant_fasta_identifier)]
    if prefix_protein == 'constant':
        chimera_sequences = [constant_sequence_of_interest[0] + variants for variants in variant_sequences]
    elif prefix_protein == 'variant':
        chimera_sequences = [variants + constant_sequence_of_interest[0] for variants in variant_sequences]
    list(map(ChimeraGenerator.fasta_creation, chimera_fastas, chimera_sequences, number_of_subunits))
    if fasta_toggles['Make a list of created fasta files?'] == '#':
        if number_of_subunits == 1:
            fasta_list = monomer_fastas + chimera_fastas + fasta_arguments['constant_fasta_for_alphafold']
        else:
            fasta_list = multimer_fastas + chimera_fastas + fasta_arguments['constant_fasta_for_alphafold']
        with open(fasta_arguments['fasta_list_file_name'], 'w') as fasta_list_file:
            for fasta in fasta_list:
                fasta_list_file.write(f'{fasta}\n')
if operation_toggles['run_analysis_operation'] == '#':
    plddt_directory = naming_arguments['plddt_directory']
    analysis_toggles = analysis_arguments['analysis_toggles']
    plddt_extension = naming_arguments['plddt_file_extension']
    naming_convention = 'monomer_naming_convention' if number_of_subunits[
                                                           0] == 1 else naming_convention = 'multimer_naming_convention'
    native_plddts = [plddt_directory + naming_arguments[naming_convention].replace(
        placeholder, protein) + plddt_extension for protein in
                     protein_names]
    chimera_plddt = [plddt_directory + naming_arguments['chimera_naming_convention'].replace(
        placeholder, protein) + plddt_extension for protein in
                     protein_names]
    plddt_files = native_plddts + chimera_plddt
    constant_plddt = [analysis_arguments['constant_plddt'] for x in chimera_plddt]
    emboss_files = [naming_arguments['emboss_names'].replace(placeholder, protein) for protein in
                    protein_names]
    pdbs = ['NA' for x in plddt_files]
    column_names = analysis_arguments['column_names']
    average_relative_stability = []
    chimera_boundary = [(0, len(constant_sequence_of_interest[0])),
                        (len(constant_sequence_of_interest[0]), None)]
    if prefix_protein == 'variant': chimera_boundary = chimera_boundary[::-1]
    alphafold_folders = [analysis_arguments['alphafold_outputs_directory'] + Path(file).stem + '/' for file in
                         plddt_files]
    if analysis_toggles['make_pdbs?'] == '#':
        pdbs = [naming_arguments['pdb_directory'] + Path(file).stem + naming_arguments['pdb_extension'] for file in
                plddt_files]
    if analysis_toggles['make_plddts?'] == '#':
        list(map(Analysis.generate_alphafold_files, alphafold_folders, plddt_files, pdbs))
        Analysis.generate_alphafold_files(
            analysis_arguments['alphafold_outputs_directory'] + analysis_arguments[
                'constant_alphafold_folder_name'] + '/',
            constant_plddt[0])
        if number_of_subunits[0] > 1:
            averaged_native_plddt = [
                plddt_directory + naming_arguments['averaged_multimer_naming_convention'].replace(placeholder,
                                                                                                  protein) + plddt_extension
                for protein
                in protein_names]
            averaged_chimera_plddt = [
                plddt_directory + naming_arguments['averaged_chimera_naming_convention'].replace(placeholder,
                                                                                                 protein) + plddt_extension
                for protein
                in protein_names]
            list(map(Analysis.averaging_multimer_plddt, native_plddts, averaged_native_plddt, number_of_subunits))
            list(map(Analysis.averaging_multimer_plddt, chimera_plddt, averaged_chimera_plddt, number_of_subunits))
            Analysis.averaging_multimer_plddt(analysis_arguments['constant_plddt'],
                                              analysis_arguments['averaged_constant'],
                                              number_of_subunits[0])
            averaged_constant_plddt = [analysis_arguments['averaged_constant'] for x in averaged_chimera_plddt]
    constant_boundary = \
    AccessiontoAlignment.alignment_finder(constant_fasta, constant_sequence_of_interest[0], constant_fasta_identifier,
                                          constant_fasta_identifier)[
        1]
    if analysis_toggles['make_emboss_files?'] == '#':
        variant_splice_info = list(
            map(AccessiontoAlignment.alignment_finder, msa, variant_sequence_of_interest, protein_names,
                variant_fasta_identifier, [analysis_arguments['emboss_command'] for x in emboss_files],
                emboss_files))
    else:
        variant_splice_info = list(
            map(AccessiontoAlignment.alignment_finder, msa, variant_sequence_of_interest, protein_names,
                variant_fasta_identifier))
    variant_splice_info = [x[1] for x in variant_splice_info]
    if number_of_subunits == 1:
        for index, boundary in enumerate(variant_splice_info):
            constant_stability = Analysis.relative_stability(constant_plddt[index], constant_boundary,
                                                             chimera_plddt[index], chimera_boundary[0])
            variant_stability = Analysis.relative_stability(native_plddts[index], variant_splice_info[index],
                                                            chimera_plddt[index], chimera_boundary[1])
            average_relative_stability.append(
                (constant_stability[0] + variant_stability[0]) / (constant_stability[1] + variant_stability[1]))
    else:
        for index, boundary in enumerate(variant_splice_info):
            constant_stability = Analysis.relative_stability(averaged_constant_plddt[index], constant_boundary,
                                                             averaged_chimera_plddt[index],
                                                             chimera_boundary[0])
            variant_stability = Analysis.relative_stability(averaged_native_plddt[index], variant_splice_info[index],
                                                            averaged_chimera_plddt[index], chimera_boundary[1])
            average_relative_stability.append(
                (constant_stability[0] + variant_stability[0]) / (constant_stability[1] + variant_stability[1]))

    if column_names['similarity'][0] == '#':
        similarity = list(map(Analysis.get_sequence_similarity, emboss_files))
    if column_names['overall_native_stability'][0] == '#':
        overall_native_stability = list(map(Analysis.overall_confidence, native_plddts))
    if column_names['overall_chimera_stability'][0] == '#':
        overall_chimera_stability = list(map(Analysis.overall_confidence, chimera_plddt))
    column_names = {key: value[1] for (key, value) in analysis_arguments['column_names'].items() if value[0] == '#'}
    data_array = empty((len(average_relative_stability) + 1, len(column_names)), dtype=object)
    for column_count, (corresponding_data, column_name) in enumerate(column_names.items()):
        exec(f'data_array[0, column_count], data_array[1:, column_count] = column_name, corresponding_data')
    savetxt(analysis_arguments['analysis_output_csv'], data_array, fmt=','.join('%s' for x in column_names),
            delimiter=",")
