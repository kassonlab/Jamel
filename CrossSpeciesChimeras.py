import ChimeraGenerator
from sys import argv
from json import load
import AccessiontoAlignment
from numpy import savetxt, empty
import Analysis
from pathlib import Path
# from setup import create_setup_slurm

argument_json = argv[1]
# TODO: Filesort and alphafold submission and gromacs simulation
#  turn list to tuples, and assign them to dictionary
data_dict={}
with open(argument_json, 'rb') as jfile:
    argument_dict = load(jfile)["arguments"]
    fasta_arguments,analysis_arguments = argument_dict['fasta_arguments'],argument_dict['analysis_arguments']
    operation_toggles = argument_dict['operation_toggles']
    prefix_protein = argument_dict['Which protein is first, constant or variant?']
    naming_arguments = argument_dict['naming_arguments']
with open(argument_dict['protein_list'], 'r') as info_list:
    info_list = info_list.readlines()
    data_dict['accession_number'],data_dict['protein_names'] = tuple(x.split()[0] for x in info_list),tuple(x.split()[1] for x in info_list)
constant_sequence_of_interest = fasta_arguments['constant_sequence_of_interest']
variant_sequence_of_interest = fasta_arguments['variant_sequence_of_interest']


with open(constant_sequence_of_interest, 'r') as fasta:
    constant_sequence_of_interest = ''.join(x for x in fasta if x[0] != '>' if x != '').strip().replace('\n', '')

with open(variant_sequence_of_interest, 'r') as fasta:
    variant_sequence_of_interest = ''.join(x for x in fasta if x[0] != '>' if x != '').strip().replace('\n', '')
list_len=len(data_dict['protein_names'])
constant_fasta,variant_fasta = fasta_arguments['constant_fasta'],fasta_arguments['variant_fasta']
placeholder = naming_arguments['placeholder']
constant_fasta_identifier = fasta_arguments['constant_fasta_identifier']
variant_fasta_identifier = fasta_arguments['variant_fasta_identifier'],
msa=fasta_arguments['msa_file_name'],
constant_sequence_of_interest = constant_sequence_of_interest
variant_sequence_of_interest = variant_sequence_of_interest,
number_of_subunits = fasta_arguments['number_of_subunits'],

if operation_toggles['run_fasta_operation'] == '#':
    fasta_data_dict=data_dict.copy()
    email = fasta_arguments['email_for_accession'],
    fasta_data_dict['monomer_fastas'] = tuple(naming_arguments['fasta_directory'] + naming_arguments['monomer_naming_convention'].replace(
        placeholder, name) + naming_arguments['fasta_file_extension'] for name in fasta_data_dict['protein_names'])
    fasta_data_dict['chimera_fastas'] = tuple(naming_arguments['fasta_directory'] + naming_arguments['chimera_naming_convention'].replace(
        placeholder, name) + naming_arguments['fasta_file_extension'] for name in fasta_data_dict['protein_names'])
    if number_of_subunits[0] == 1:
        tuple(map(AccessiontoAlignment.accession_to_fasta, fasta_data_dict['monomer_fastas'], fasta_data_dict['accession_number'], email*list_len, number_of_subunits*list_len))
    else:
        fasta_data_dict['multimer_fastas'] = tuple(naming_arguments['fasta_directory'] + naming_arguments['multimer_naming_convention'].replace(
            placeholder, name) + naming_arguments['fasta_file_extension'] for name in fasta_data_dict['protein_names'])
        tuple(map(AccessiontoAlignment.accession_to_fasta, fasta_data_dict['monomer_fastas'], fasta_data_dict['accession_number'], email*list_len, number_of_subunits*list_len,
                 fasta_data_dict['multimer_fastas']))
    fasta_toggles = fasta_arguments['fasta_toggles']
    if fasta_toggles['Create an alignment?'] == '#':
        AccessiontoAlignment.multiple_sequence_alignment(fasta_data_dict['monomer_fastas'], fasta_arguments['msa_fasta'], msa[0],
                                                         variant_fasta, fasta_arguments['muscle_command_for_msa'])
    fasta_data_dict['variant_sequences'] = tuple(x[0] for x in
                         map(AccessiontoAlignment.alignment_finder, msa*list_len, variant_sequence_of_interest*list_len, fasta_data_dict['protein_names'],
                             variant_fasta_identifier*list_len))
    if prefix_protein == 'constant':
        fasta_data_dict['chimera_sequences'] = tuple(constant_sequence_of_interest + variants for variants in fasta_data_dict['variant_sequences'])
    elif prefix_protein == 'variant':
        fasta_data_dict['chimera_sequences'] = tuple(variants + constant_sequence_of_interest for variants in fasta_data_dict['variant_sequences'])
    tuple(map(ChimeraGenerator.fasta_creation, fasta_data_dict['chimera_fastas'], fasta_data_dict['chimera_sequences'], number_of_subunits*list_len))
    if fasta_toggles['Make a list of created fasta files?'] == '#':
        if number_of_subunits[0] == 1:
            fasta_data_dict['fasta_list'] = tuple(fasta_data_dict['monomer_fastas']) + tuple(fasta_data_dict['chimera_fastas']) + (fasta_arguments['constant_fasta_for_alphafold'],)
        else:
            fasta_data_dict['fasta_list'] = tuple(fasta_data_dict['multimer_fastas']) + tuple(fasta_data_dict['chimera_fastas']) + (fasta_arguments['constant_fasta_for_alphafold'],)
        with open(fasta_arguments['fasta_list_file_name'], 'w') as fasta_list_file:
            for fasta in fasta_data_dict['fasta_list']:
                fasta_list_file.write(f'{fasta}\n')
    del fasta_data_dict
if operation_toggles['run_analysis_operation'] == '#':
    analysis_data_dict=data_dict.copy()
    plddt_directory,plddt_extension = naming_arguments['plddt_directory'],naming_arguments['plddt_file_extension']
    analysis_toggles = analysis_arguments['analysis_toggles']
    naming_convention = 'monomer_naming_convention' if number_of_subunits[
                                                           0] == 1 else 'multimer_naming_convention'
    analysis_data_dict['native_plddts'] = tuple(plddt_directory + naming_arguments[naming_convention].replace(
        placeholder, protein) + plddt_extension for protein in
                     analysis_data_dict['protein_names'])
    analysis_data_dict['chimera_plddts'] = tuple(plddt_directory + naming_arguments['chimera_naming_convention'].replace(
        placeholder, protein) + plddt_extension for protein in
                           analysis_data_dict['protein_names'])
    analysis_data_dict['plddt_files'] = analysis_data_dict['native_plddts'] + analysis_data_dict['chimera_plddts']
    constant_plddt = analysis_arguments['constant_plddt']
    emboss_files = tuple(naming_arguments['emboss_names'].replace(placeholder, protein) for protein in
                    analysis_data_dict['protein_names'])
    chimera_boundary = [(0, len(constant_sequence_of_interest)),
                        (len(constant_sequence_of_interest), None)]
    if prefix_protein == 'variant': chimera_boundary = chimera_boundary[::-1]
    analysis_data_dict['alphafold_folders'] = tuple(analysis_arguments['alphafold_outputs_directory'] + Path(file).stem + '/' for file in
                         analysis_data_dict['plddt_files'])

    analysis_data_dict['pdbs'] = tuple('NA' for x in analysis_data_dict['plddt_files'])
    if analysis_toggles['make_pdbs?'] == '#':
        analysis_data_dict['pdbs'] = tuple(naming_arguments['pdb_directory'] + Path(file).stem + naming_arguments['pdb_extension'] for file in
                analysis_data_dict['plddt_files'])
    if number_of_subunits[0]>1:
        analysis_data_dict['averaged_native_plddt'] = tuple(
            plddt_directory + naming_arguments['averaged_multimer_naming_convention'].replace(placeholder,
                                                                                              protein) + plddt_extension
            for protein
            in analysis_data_dict['protein_names'])
        analysis_data_dict['averaged_chimera_plddt'] = tuple(
            plddt_directory + naming_arguments['averaged_chimera_naming_convention'].replace(placeholder,
                                                                                             protein) + plddt_extension
            for protein
            in analysis_data_dict['protein_names'])
        averaged_constant_plddt = analysis_arguments['averaged_constant']
    if analysis_toggles['make_plddts?'] == '#':
        tuple(map(Analysis.generate_alphafold_files, analysis_data_dict['alphafold_folders'], analysis_data_dict['plddt_files'], analysis_data_dict['pdbs']))
        Analysis.generate_alphafold_files(
            analysis_arguments['alphafold_outputs_directory'] + analysis_arguments[
                'constant_alphafold_folder_name'] + '/',
            constant_plddt)
        if number_of_subunits[0] > 1:
            tuple(map(Analysis.averaging_multimer_plddt, analysis_data_dict['native_plddts'], analysis_data_dict['averaged_native_plddt'], number_of_subunits*list_len))
            tuple(map(Analysis.averaging_multimer_plddt, analysis_data_dict['chimera_plddts'], analysis_data_dict['averaged_chimera_plddt'], number_of_subunits * list_len))
            Analysis.averaging_multimer_plddt(constant_plddt,
                                              analysis_arguments['averaged_constant'],
                                              number_of_subunits[0])
    if analysis_toggles['make_emboss_files?'] == '#':
        variant_splice_info = tuple(
            map(AccessiontoAlignment.alignment_finder, msa*list_len, variant_sequence_of_interest*list_len, analysis_data_dict['protein_names'],
                variant_fasta_identifier*list_len, [analysis_arguments['emboss_command'] for x in emboss_files],
                emboss_files))
    else:
        variant_splice_info = tuple(
            map(AccessiontoAlignment.alignment_finder, msa*list_len, variant_sequence_of_interest*list_len, analysis_data_dict['protein_names'],
                variant_fasta_identifier*list_len))
    analysis_data_dict['variant_splice_info'] = tuple(boundary[1] for boundary in variant_splice_info)
    constant_boundary = \
        AccessiontoAlignment.alignment_finder(constant_fasta, constant_sequence_of_interest, constant_fasta_identifier,
                                              constant_fasta_identifier)[
            1]
    analysis_data_dict['average_relative_stability'] = []
    if number_of_subunits[0] == 1:
        for index, boundary in enumerate(analysis_data_dict['variant_splice_info']):
            constant_stability = Analysis.relative_stability(constant_plddt, constant_boundary,
                                                             analysis_data_dict['chimera_plddts'][index], chimera_boundary[0])
            variant_stability = Analysis.relative_stability(analysis_data_dict['native_plddts'][index], analysis_data_dict['variant_splice_info'][index],
                                                            analysis_data_dict['chimera_plddts'][index], chimera_boundary[1])
            analysis_data_dict['average_relative_stability'].append(
                (constant_stability[0] + variant_stability[0]) / (constant_stability[1] + variant_stability[1]))
    else:
        for index, boundary in enumerate(analysis_data_dict['variant_splice_info']):
            constant_stability = Analysis.relative_stability(averaged_constant_plddt, constant_boundary,
                                                             analysis_data_dict['averaged_chimera_plddt'][index],
                                                             chimera_boundary[0])
            variant_stability = Analysis.relative_stability(analysis_data_dict['averaged_native_plddt'][index], analysis_data_dict['variant_splice_info'][index],
                                                            analysis_data_dict['averaged_chimera_plddt'][index], chimera_boundary[1])
            analysis_data_dict['average_relative_stability'].append(
                (constant_stability[0] + variant_stability[0]) / (constant_stability[1] + variant_stability[1]))
    column_names = analysis_arguments['column_names']
    if column_names['similarity'][0] == '#' and analysis_toggles['make_emboss_files?']=="#":
        analysis_data_dict['similarity'] = tuple(map(Analysis.get_sequence_similarity, emboss_files))
    if column_names['overall_native_stability'][0] == '#':
        analysis_data_dict['overall_native_stability'] = tuple(map(Analysis.overall_confidence, analysis_data_dict['native_plddts']))
    if column_names['overall_chimera_stability'][0] == '#':
        analysis_data_dict['overall_chimera_stability'] = tuple(map(Analysis.overall_confidence, analysis_data_dict['chimera_plddts']))
    column_names = {key: value[1] for (key, value) in analysis_arguments['column_names'].items() if value[0] == '#'}
    data_array = empty((len(analysis_data_dict['average_relative_stability']) + 1, len(column_names)), dtype=object)
    for column_count, (corresponding_data, column_name) in enumerate(column_names.items()):
        exec(f'data_array[0, column_count], data_array[1:, column_count] = column_name, corresponding_data')
    savetxt(analysis_arguments['analysis_output_csv'], data_array, fmt=','.join('%s' for x in column_names),
            delimiter=",")
    del analysis_data_dict
    # TODO have a functiuon to proint pdbs
# TODO: Create a way for python to check periodically whether the setup slurms are done and then submit the prodcution
if operation_toggles['run_gromacs_operation']=='#':
    gromacs_arguments = argument_dict['gromacs_arguments']
    gromacs_toggles = gromacs_arguments['gromacs_toggles']
    with open(gromacs_arguments['pdbs_to_run'],'r') as pdbs_to_run:
        pdbs_to_run=(pdb.split()[0] for pdb in pdbs_to_run.readlines())
    if gromacs_toggles['create setup slurms']=='#':
        ''