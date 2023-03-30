import ChimeraGenerator
from sys import argv
from json import load
import AccessiontoAlignment
from numpy import savetxt, empty
import Analysis
from pathlib import Path
from os import system,path
from csv import reader
from setup import create_setup_slurm,create_prod_slurm

argument_json = argv[1]
with open(argument_json, 'rb') as jfile:
    argument_dict = load(jfile)["arguments"]
    fasta_arguments,analysis_arguments = argument_dict['fasta_arguments'],argument_dict['analysis_arguments']
    operation_toggles = argument_dict['operation_toggles']
    prefix_protein = argument_dict['Which protein is first, constant or variant?']
    naming_arguments = argument_dict['naming_arguments']
chimera_container=()
with open(argument_dict['protein_list'], 'r') as info_list:
    info_list = info_list.readlines()
    for index,line in enumerate(info_list):
        exec(f'chimera_{index}=ChimeraGenerator.chimeracls()')
        chimera_container += globals()[f'chimera_{index}'],
        globals()[f'chimera_{index}'].nickname=line.split()[1]
        globals()[f'chimera_{index}'].accession = line.split()[0]
constant_sequence_of_interest = fasta_arguments['constant_sequence_of_interest']
variant_sequence_of_interest = fasta_arguments['variant_sequence_of_interest']
with open(constant_sequence_of_interest, 'r') as fasta:
    constant_sequence_of_interest = ''.join(x for x in fasta if x[0] != '>' if x != '').strip().replace('\n', '')
with open(variant_sequence_of_interest, 'r') as fasta:
    variant_sequence_of_interest = ''.join(x for x in fasta if x[0] != '>' if x != '').strip().replace('\n', '')
num_of_chi=len(chimera_container)
constant_fasta,variant_fasta = fasta_arguments['constant_fasta'],fasta_arguments['variant_fasta']
placeholder = naming_arguments['placeholder']
constant_fasta_identifier = fasta_arguments['constant_fasta_identifier']
variant_fasta_identifier = fasta_arguments['variant_fasta_identifier']
msa=fasta_arguments['msa_file_name']
subunits = fasta_arguments['number_of_subunits']
for chimera in chimera_container:
    chimera.monomer_stem=naming_arguments['monomer_naming_convention'].replace(placeholder, chimera.nickname)
    chimera.chimera_stem=naming_arguments['chimera_naming_convention'].replace(placeholder, chimera.nickname)
    chimera.chi_pdb = naming_arguments['pdb_directory'] + chimera.chimera_stem + naming_arguments['pdb_extension']
    if subunits>1:
        chimera.multimer_stem=naming_arguments['multimer_naming_convention'].replace(placeholder, chimera.nickname)
        chimera.native_pdb=naming_arguments['pdb_directory'] + chimera.multimer_stem + naming_arguments['pdb_extension']
    else:
        chimera.native_pdb = naming_arguments['pdb_directory'] + chimera.monomer_stem + naming_arguments['pdb_extension']
if operation_toggles['run_fasta_operation'] == '#':
    email = fasta_arguments['email_for_accession']
    fasta_toggles = fasta_arguments['fasta_toggles']
    for chimera in chimera_container:
        monomer_fasta = naming_arguments['fasta_directory'] + chimera.monomer_stem + naming_arguments['fasta_file_extension']
        chimera_fastas = naming_arguments['fasta_directory'] + chimera.chimera_stem + naming_arguments['fasta_file_extension']
        if subunits == 1:
            AccessiontoAlignment.accession_to_fasta(monomer_fasta, chimera.accession, (email,) * num_of_chi, (subunits,) * num_of_chi)
        else:
            multimer_fasta=naming_arguments['fasta_directory'] + chimera.multimer_stem + naming_arguments['fasta_file_extension']
            AccessiontoAlignment.accession_to_fasta(monomer_fasta, chimera.accession, (email,) * num_of_chi, (subunits,) * num_of_chi,
                      multimer_fasta)
    if fasta_toggles['Create an alignment?'] == '#':
        AccessiontoAlignment.multiple_sequence_alignment(tuple(naming_arguments['fasta_directory'] + chimera.monomer_stem + naming_arguments['fasta_file_extension'] for chimera in chimera_container), fasta_arguments['msa_fasta'], msa,
                                                         variant_fasta, fasta_arguments['muscle_command_for_msa'])
    for chimera in chimera_container:
        variant_sequence=AccessiontoAlignment.alignment_finder(msa, variant_sequence_of_interest, chimera.nickname,
                             variant_fasta_identifier)[0]
        if prefix_protein == 'constant':
            chimera_sequence = constant_sequence_of_interest + variant_sequence
        elif prefix_protein == 'variant':
            chimera_sequence = variant_sequence + constant_sequence_of_interest
        chimera_fasta=naming_arguments['fasta_directory'] + chimera.chimera_stem + naming_arguments['fasta_file_extension']
        ChimeraGenerator.fasta_creation(chimera_fasta, chimera_sequence, subunits)
        if fasta_toggles['Make a list of created fasta files?'] == '#':
            with open(fasta_arguments['fasta_list_file_name'], 'w') as fasta_list_file:
                fasta_list_file.write(f'{chimera_fasta}\n')
                if subunits == 1:
                    fasta_list_file.write(f"{naming_arguments['fasta_directory'] + chimera.monomer_stem + naming_arguments['fasta_file_extension']}\n")
                else:
                    fasta_list_file.write(
                        f"{naming_arguments['fasta_directory'] + chimera.multimer_stem + naming_arguments['fasta_file_extension']}\n")
    if fasta_toggles['Make a list of created fasta files?'] == '#':
        with open(fasta_arguments['fasta_list_file_name'], 'w') as fasta_list_file:
            fasta_list_file.write(fasta_arguments['constant_fasta_for_alphafold'])
if operation_toggles['run_analysis_operation'] == '#':
    plddt_directory,plddt_extension = naming_arguments['plddt_directory'],naming_arguments['plddt_file_extension']
    analysis_toggles = analysis_arguments['analysis_toggles']
    constant_plddt = analysis_arguments['constant_plddt']
    chimera_boundary = [(0, len(constant_sequence_of_interest)),
                        (len(constant_sequence_of_interest), None)]
    if prefix_protein == 'variant': chimera_boundary = chimera_boundary[::-1]
    averaged_constant_plddt = analysis_arguments['averaged_constant']
    constant_boundary = \
        AccessiontoAlignment.alignment_finder(constant_fasta, constant_sequence_of_interest, constant_fasta_identifier,
                                              constant_fasta_identifier)[
            1]
    if analysis_toggles['make_plddts?'] == '#':
        Analysis.generate_alphafold_files(
            analysis_arguments['alphafold_outputs_directory'] + analysis_arguments[
                'constant_alphafold_folder_name'] + '/',
            constant_plddt)
        if subunits > 1:
            Analysis.averaging_multimer_plddt(constant_plddt,
                                              analysis_arguments['averaged_constant'],
                                              subunits)
    for chimera in chimera_container:
        chimera.native_plddt=plddt_directory + chimera.monomer_stem + plddt_extension if subunits==1 else plddt_directory + chimera.multimer_stem + plddt_extension
        chimera.chi_plddt=plddt_directory + chimera.chimera_stem + plddt_extension
        native_folder=analysis_arguments['alphafold_outputs_directory'] + chimera.monomer_stem + '/' if subunits==1 else analysis_arguments['alphafold_outputs_directory'] + chimera.multimer_stem + '/'
        chimera_folder=analysis_arguments['alphafold_outputs_directory'] + chimera.chimera_stem + '/'
        if analysis_toggles['make_pdbs?'] == '#':
            Analysis.generate_alphafold_files(native_folder,"NA",chimera.native_pdb)
            Analysis.generate_alphafold_files(chimera_folder, "NA", chimera.chi_pdb)
        if subunits>1:
            avg_native_plddt = plddt_directory + naming_arguments['averaged_multimer_naming_convention'].replace(placeholder, chimera.nickname) + plddt_extension
            avg_chimera_plddt = plddt_directory + naming_arguments['averaged_chimera_naming_convention'].replace(placeholder,
                                                                                                 chimera.nickname) + plddt_extension
        if analysis_toggles['make_plddts?'] == '#':
            Analysis.generate_alphafold_files(native_folder, chimera.native_plddt)
            Analysis.generate_alphafold_files(chimera_folder, chimera.chi_plddt)
            if subunits > 1:
                Analysis.averaging_multimer_plddt(chimera.native_plddt, avg_native_plddt, subunits )
                Analysis.averaging_multimer_plddt(chimera.chi_plddt, avg_chimera_plddt, subunits)
        if analysis_toggles['make_emboss_files?'] == '#':
            emboss=naming_arguments['emboss_names'].replace(placeholder, chimera.nickname)
            variant_splice_info =AccessiontoAlignment.alignment_finder(msa, variant_sequence_of_interest, chimera.nickname,
                    variant_fasta_identifier, analysis_arguments['emboss_command'],
                    emboss)[1]
        else:
            variant_splice_info = AccessiontoAlignment.alignment_finder(msa, variant_sequence_of_interest, chimera.nickname,
                    variant_fasta_identifier)[1]
        if subunits == 1:
            constant_stability = Analysis.relative_stability(constant_plddt, constant_boundary,
                                                             chimera.chi_plddt, chimera_boundary[0])
            variant_stability = Analysis.relative_stability(chimera.native_plddt, variant_splice_info,
                                                            chimera.chi_plddt, chimera_boundary[1])
            chimera.rel_stability=(constant_stability[0] + variant_stability[0]) / (constant_stability[1] + variant_stability[1])
        else:
            constant_stability = Analysis.relative_stability(averaged_constant_plddt, constant_boundary,
                                                             avg_chimera_plddt,
                                                             chimera_boundary[0])
            variant_stability = Analysis.relative_stability(avg_native_plddt, variant_splice_info,
                                                            avg_chimera_plddt, chimera_boundary[1])
            chimera.rel_stability=(constant_stability[0] + variant_stability[0]) / (constant_stability[1] + variant_stability[1])
    column_names = analysis_arguments['column_names']
    data_dict= {'average_relative_stability': tuple(chimera.rel_stability for chimera in chimera_container),'nickname': tuple(chimera.nickname for chimera in chimera_container)}
    if column_names['similarity'][0] == '#':
        emboss_files = tuple(naming_arguments['emboss_names'].replace(placeholder, chimera.nickname) for chimera in
                             chimera_container)
        data_dict['similarity'] = tuple(map(Analysis.get_sequence_similarity, emboss_files))
    if column_names['overall_native_stability'][0] == '#':
        data_dict['overall_native_stability'] = tuple(Analysis.overall_confidence(chimera.native_plddt) for chimera in chimera_container)
    if column_names['overall_chimera_stability'][0] == '#':
        data_dict['overall_chimera_stability'] = tuple(Analysis.overall_confidence(chimera.chi_plddt) for chimera in chimera_container)
    column_names = {key: value[1] for (key, value) in analysis_arguments['column_names'].items() if value[0] == '#'}
    data_array = empty((len(data_dict['average_relative_stability']) + 1, len(column_names)), dtype=object)
    for column_count, (corresponding_data, column_name) in enumerate(column_names.items()):
        data_array[0, column_count], data_array[1:, column_count] = column_name, data_dict[corresponding_data]
    savetxt(analysis_arguments['analysis_output_csv'], data_array, fmt=','.join('%s' for x in column_names),
            delimiter=",")
    del data_dict
    # TODO have a functiuon to proint pdbs and make pdbs for reference
    # TODO make it to be able to simulate native pdb
# TODO: Create a way for python to check periodically whether the setup slurms are done and then submit the prodcution
# TODO make certain parts available from commandline
if operation_toggles['run_gromacs_operation']=='#':
    gromacs_data_dict={}
    gromacs_arguments = argument_dict['gromacs_arguments']
    gromacs_toggles = gromacs_arguments['gromacs_toggles']
    # TODO if you hyand select make a file to resubmit with the choisces
    if not path.exists(gromacs_arguments['pdbs_to_run']) or gromacs_toggles['create new pdb list']:
        with open(analysis_arguments['analysis_output_csv']) as data:
            pdbs_to_run=tuple(chimera.chi_pdb for chimera in chimera_container)
            for row in reader(data):
                print(row)
            for index,pdb in enumerate(pdbs_to_run):
                print(f'{index}. {pdb}')
            pdbs_to_run=tuple(pdbs_to_run[int(index)] for index in input("Submit comma separated list:\n Ex: 0,1,2,7,8\n Enter:").split(','))
            with open(gromacs_arguments['pdbs_to_run'],'w') as pdb_list:
                for pdb in pdbs_to_run:
                    pdb_list.write(f'{pdb}\n')
    else:
        with open(gromacs_arguments['pdbs_to_run'],'r') as pdbs_to_run:
            pdbs_to_run=tuple(pdb.split()[0] for pdb in pdbs_to_run.readlines())
    for pdb in pdbs_to_run:
        system(f'cp {pdb} {naming_arguments["gromacs_slurm_dir"]}')
    gromacs_data_dict['setup_slurms'] = tuple(
        naming_arguments['gromacs_slurm_dir'] + Path(pdb).stem + naming_arguments['setup_extension'] for pdb in
        pdbs_to_run)
    # TODO make folders for gromacs
    if gromacs_toggles['create setup slurms']=='#':
        for index,(pdb,slurm) in enumerate(zip(pdbs_to_run,gromacs_data_dict['setup_slurms'])):
            create_setup_slurm(pdb,gromacs_arguments['gmxbin'],gromacs_arguments['pdb2gmx'],gromacs_arguments['slurm_template'],slurm,gromacs_arguments['slurm_output'].replace(placeholder,Path(pdb).stem),gromacs_arguments['slurm_error'].replace(placeholder,Path(pdb).stem))
    if gromacs_toggles['sbatch setup slurms'] == '#':
        for slurm in gromacs_data_dict['setup_slurms']:
            system(f"sbatch {slurm}")
    gromacs_data_dict['mdrun_slurms'] = tuple(
        naming_arguments['gromacs_slurm_dir'] + Path(pdb).stem + naming_arguments['production_extension'] for pdb
        in pdbs_to_run)
    if gromacs_toggles['create mdrun slurms'] == '#':
        for index, (pdb, slurm) in enumerate(zip(pdbs_to_run, gromacs_data_dict['mdrun_slurms'])):
            create_prod_slurm(pdb,gromacs_arguments['gmxbin'],gromacs_arguments['slurm_template'],slurm,gromacs_arguments['slurm_output'].replace(placeholder,Path(pdb).stem),gromacs_arguments['slurm_error'].replace(placeholder,Path(pdb).stem))
    if gromacs_toggles['sbatch mdrun slurms'] == '#':
        for slurm in gromacs_data_dict['mdrun_slurms']:
            system(f"sbatch {slurm}")
