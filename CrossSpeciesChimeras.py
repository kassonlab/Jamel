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
#  turn list to tuples
with open(argument_json, 'rb') as jfile:
    argument_dict = load(jfile)["arguments"]
    fasta_arguments,analysis_arguments = argument_dict['fasta_arguments'],argument_dict['analysis_arguments']
    operation_toggles = argument_dict['operation_toggles']
    prefix_protein = argument_dict['Which protein is first, constant or variant?']
    naming_arguments = argument_dict['naming_arguments']
with open(argument_dict['protein_list'], 'r') as info_list:
    info_list = info_list.readlines()
    accession_number,protein_names = tuple(x.split()[0] for x in info_list),tuple(x.split()[1] for x in info_list)
constant_sequence_of_interest = fasta_arguments['constant_sequence_of_interest']
variant_sequence_of_interest = fasta_arguments['variant_sequence_of_interest']


with open(constant_sequence_of_interest, 'r') as fasta:
    constant_sequence_of_interest = ''.join(x for x in fasta if x[0] != '>' if x != '').strip().replace('\n', '')

with open(variant_sequence_of_interest, 'r') as fasta:
    variant_sequence_of_interest = ''.join(x for x in fasta if x[0] != '>' if x != '').strip().replace('\n', '')
list_len=len(protein_names)
constant_fasta,variant_fasta = fasta_arguments['constant_fasta'],fasta_arguments['variant_fasta']
placeholder = naming_arguments['placeholder']
constant_fasta_identifier = fasta_arguments['constant_fasta_identifier']
variant_fasta_identifier = fasta_arguments['variant_fasta_identifier'],
msa=fasta_arguments['msa_file_name'],
constant_sequence_of_interest = constant_sequence_of_interest
variant_sequence_of_interest = variant_sequence_of_interest,
number_of_subunits = fasta_arguments['number_of_subunits'],

if operation_toggles['run_fasta_operation'] == '#':
    email = fasta_arguments['email_for_accession'],
    monomer_fastas = tuple(naming_arguments['fasta_directory'] + naming_arguments['monomer_naming_convention'].replace(
        placeholder, name) + naming_arguments['fasta_file_extension'] for name in protein_names)
    chimera_fastas = tuple(naming_arguments['fasta_directory'] + naming_arguments['chimera_naming_convention'].replace(
        placeholder, name) + naming_arguments['fasta_file_extension'] for name in protein_names)
    if number_of_subunits[0] == 1:
        tuple(map(AccessiontoAlignment.accession_to_fasta, monomer_fastas, accession_number, email*list_len, number_of_subunits*list_len))
    else:
        multimer_fastas = tuple(naming_arguments['fasta_directory'] + naming_arguments['multimer_naming_convention'].replace(
            placeholder, name) + naming_arguments['fasta_file_extension'] for name in protein_names)
        tuple(map(AccessiontoAlignment.accession_to_fasta, monomer_fastas, accession_number, email*list_len, number_of_subunits*list_len,
                 multimer_fastas))
    fasta_toggles = fasta_arguments['fasta_toggles']
    if fasta_toggles['Create an alignment?'] == '#':
        AccessiontoAlignment.multiple_sequence_alignment(monomer_fastas, fasta_arguments['msa_fasta'], msa[0],
                                                         variant_fasta, fasta_arguments['muscle_command_for_msa'])
    variant_sequences = tuple(x[0] for x in
                         map(AccessiontoAlignment.alignment_finder, msa*list_len, variant_sequence_of_interest*list_len, protein_names,
                             variant_fasta_identifier*list_len))
    if prefix_protein == 'constant':
        chimera_sequences = tuple(constant_sequence_of_interest + variants for variants in variant_sequences)
    elif prefix_protein == 'variant':
        chimera_sequences = tuple(variants + constant_sequence_of_interest for variants in variant_sequences)
    tuple(map(ChimeraGenerator.fasta_creation, chimera_fastas, chimera_sequences, number_of_subunits*list_len))
    if fasta_toggles['Make a list of created fasta files?'] == '#':
        if number_of_subunits[0] == 1:
            fasta_list = tuple(monomer_fastas) + tuple(chimera_fastas) + (fasta_arguments['constant_fasta_for_alphafold'],)
        else:
            fasta_list = tuple(multimer_fastas) + tuple(chimera_fastas) + (fasta_arguments['constant_fasta_for_alphafold'],)
        with open(fasta_arguments['fasta_list_file_name'], 'w') as fasta_list_file:
            for fasta in fasta_list:
                fasta_list_file.write(f'{fasta}\n')
if operation_toggles['run_analysis_operation'] == '#':
    plddt_directory,plddt_extension = naming_arguments['plddt_directory'],naming_arguments['plddt_file_extension']
    analysis_toggles = analysis_arguments['analysis_toggles']
    naming_convention = 'monomer_naming_convention' if number_of_subunits[
                                                           0] == 1 else 'multimer_naming_convention'
    native_plddts = tuple(plddt_directory + naming_arguments[naming_convention].replace(
        placeholder, protein) + plddt_extension for protein in
                     protein_names)
    chimera_plddt = tuple(plddt_directory + naming_arguments['chimera_naming_convention'].replace(
        placeholder, protein) + plddt_extension for protein in
                     protein_names)
    plddt_files = native_plddts + chimera_plddt
    constant_plddt = analysis_arguments['constant_plddt']
    emboss_files = tuple(naming_arguments['emboss_names'].replace(placeholder, protein) for protein in
                    protein_names)
    chimera_boundary = [(0, len(constant_sequence_of_interest)),
                        (len(constant_sequence_of_interest), None)]
    if prefix_protein == 'variant': chimera_boundary = chimera_boundary[::-1]
    alphafold_folders = tuple(analysis_arguments['alphafold_outputs_directory'] + Path(file).stem + '/' for file in
                         plddt_files)

    pdbs = tuple('NA' for x in plddt_files)
    if analysis_toggles['make_pdbs?'] == '#':
        pdbs = tuple(naming_arguments['pdb_directory'] + Path(file).stem + naming_arguments['pdb_extension'] for file in
                plddt_files)
    if number_of_subunits[0]>1:
        averaged_native_plddt = tuple(
            plddt_directory + naming_arguments['averaged_multimer_naming_convention'].replace(placeholder,
                                                                                              protein) + plddt_extension
            for protein
            in protein_names)
        averaged_chimera_plddt = tuple(
            plddt_directory + naming_arguments['averaged_chimera_naming_convention'].replace(placeholder,
                                                                                             protein) + plddt_extension
            for protein
            in protein_names)
        print(averaged_chimera_plddt)
        averaged_constant_plddt = analysis_arguments['averaged_constant']
    if analysis_toggles['make_plddts?'] == '#':
        tuple(map(Analysis.generate_alphafold_files, alphafold_folders, plddt_files, pdbs))
        Analysis.generate_alphafold_files(
            analysis_arguments['alphafold_outputs_directory'] + analysis_arguments[
                'constant_alphafold_folder_name'] + '/',
            constant_plddt)
        if number_of_subunits[0] > 1:

            tuple(map(Analysis.averaging_multimer_plddt, native_plddts, averaged_native_plddt, number_of_subunits*list_len))
            tuple(map(Analysis.averaging_multimer_plddt, chimera_plddt, averaged_chimera_plddt, number_of_subunits*list_len))
            Analysis.averaging_multimer_plddt(constant_plddt,
                                              analysis_arguments['averaged_constant'],
                                              number_of_subunits[0])
    if analysis_toggles['make_emboss_files?'] == '#':
        variant_splice_info = tuple(
            map(AccessiontoAlignment.alignment_finder, msa*list_len, variant_sequence_of_interest*list_len, protein_names,
                variant_fasta_identifier*list_len, [analysis_arguments['emboss_command'] for x in emboss_files],
                emboss_files))
    else:
        variant_splice_info = tuple(
            map(AccessiontoAlignment.alignment_finder, msa*list_len, variant_sequence_of_interest*list_len, protein_names,
                variant_fasta_identifier*list_len))
    variant_splice_info = tuple(boundary[1] for boundary in variant_splice_info)
    constant_boundary = \
        AccessiontoAlignment.alignment_finder(constant_fasta, constant_sequence_of_interest, constant_fasta_identifier,
                                              constant_fasta_identifier)[
            1]
    average_relative_stability = []
    if number_of_subunits[0] == 1:
        for index, boundary in enumerate(variant_splice_info):
            constant_stability = Analysis.relative_stability(constant_plddt, constant_boundary,
                                                             chimera_plddt[index], chimera_boundary[0])
            variant_stability = Analysis.relative_stability(native_plddts[index], variant_splice_info[index],
                                                            chimera_plddt[index], chimera_boundary[1])
            average_relative_stability.append(
                (constant_stability[0] + variant_stability[0]) / (constant_stability[1] + variant_stability[1]))
    else:
        for index, boundary in enumerate(variant_splice_info):
            constant_stability = Analysis.relative_stability(averaged_constant_plddt, constant_boundary,
                                                             averaged_chimera_plddt[index],
                                                             chimera_boundary[0])
            variant_stability = Analysis.relative_stability(averaged_native_plddt[index], variant_splice_info[index],
                                                            averaged_chimera_plddt[index], chimera_boundary[1])
            average_relative_stability.append(
                (constant_stability[0] + variant_stability[0]) / (constant_stability[1] + variant_stability[1]))
    column_names = analysis_arguments['column_names']
    if column_names['similarity'][0] == '#' and analysis_toggles['make_emboss_files?']=="#":
        similarity = tuple(map(Analysis.get_sequence_similarity, emboss_files))
    if column_names['overall_native_stability'][0] == '#':
        overall_native_stability = tuple(map(Analysis.overall_confidence, native_plddts))
    if column_names['overall_chimera_stability'][0] == '#':
        overall_chimera_stability = tuple(map(Analysis.overall_confidence, chimera_plddt))
    column_names = {key: value[1] for (key, value) in analysis_arguments['column_names'].items() if value[0] == '#'}
    data_array = empty((len(average_relative_stability) + 1, len(column_names)), dtype=object)
    for column_count, (corresponding_data, column_name) in enumerate(column_names.items()):
        exec(f'data_array[0, column_count], data_array[1:, column_count] = column_name, corresponding_data')
    savetxt(analysis_arguments['analysis_output_csv'], data_array, fmt=','.join('%s' for x in column_names),
            delimiter=",")
    # TODO have a functiuon to proint pdbs
# TODO: Create a way for python to check periodically whether the setup slurms are done and then submit the prodcution
if operation_toggles['run_gromacs_operation']=='#':
    gromacs_arguments = argument_dict['gromacs_arguments']
    gromacs_toggles = gromacs_arguments['gromacs_toggles']
    with open(gromacs_arguments['pdbs_to_run'],'r') as pdbs_to_run:
        pdbs_to_run=(pdb.split()[0] for pdb in pdbs_to_run.readlines())
    if gromacs_toggles['create setup slurms']=='#':
        ''