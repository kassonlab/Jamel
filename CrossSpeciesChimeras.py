import copy
from json import load
from os import system, path
from pathlib import Path
from sys import exit
from time import perf_counter
from numpy import empty, savetxt
import Analysis
from setup import create_alphafold_slurm
from AccessiontoAlignment import alignment_finder,accession_to_fasta,multiple_sequence_alignment
from ChimeraGenerator import fasta_creation, chimeracls, update_json
from itertools import product
import argparse
# TODO add a commandline to change specific keys and replace recurring string through out all keys
# and be able to request operation from command line that overwrites whats in json
parser = argparse.ArgumentParser(
    description='Creating chimeric proteins from one or more json proteins corresponding to each unique protein in '
                'the structure')
parser.add_argument('-u', '--updatejson', type=str, required=False,
                    help='updating json configs ex. default_json,old_json')
parser.add_argument('-i', '--jsoninput', dest='arg_jsons', required=False, type=str,
                    help='Json config file input, comma separated for each unique protein. Make sure initial json has '
                         'naming conventions')
args = parser.parse_args()

if args.updatejson:
    args.updatejson = args.updatejson.split(',')
    update_json(args.updatejson[0], args.updatejson[1])
    exit()

class CrossNamingArguments:
    monomer_naming_convention: str
    '''A combination of words and placeholders to create a naming that will persist through all monomer chimera related data files e.g. (fasta,pdb,plddt)'''
    multimer_naming_convention: str
    chimera_naming_convention: str
    gromacs_slurm_dir: str
    emboss_names: str
    setup_extension: str
    production_extension: str
    placeholder: str
    '''a unique character to be replaced by the reference fasta identifier'''
    fasta_directory: str
    plddt_directory: str
    pdb_directory: str
    alphafold_outputs_dir: str
    fasta_extension: str
    plddt_extension: str
    pdb_extension: str
    slurm_naming: str

    def __init__(self, naming_dict):
        for key, value in naming_dict.items():
            setattr(self, key, value)


class CrossFastaArguments:
    fasta_toggles: dict
    '''A dictionary that holds optional operations within the fasta operation, they are generally turned on by 
    placing True within the quotes of their value'''
    Make_a_lst_of_created_fasta_files: str
    '''Create a file with line separted list of all fastas created'''
    Create_an_alignment: str
    constant_or_variant: str
    '''Distinguishes whether the config json is for the constant part of the chimera or the variant protein that will be swapped with all protein in the protein list'''

    protein_list: str
    '''Two column list of accession number and their associated protein nicknames'''
    fasta_identifier: str
    '''Fasta label for the constant protein or the reference protein for the homologs to be iterated through'''
    # TODO add defaults like this
    muscle_command_for_msa: str = 'module load gcc/9.2.0 && module load muscle/3.8.31 && muscle'
    '''The start of the muscle command to be run to create the alignment.'''
    number_of_subunits: int
    sequence_of_interest: str
    '''This is the path to the fasta containing the sequence of either the constant protein you want spliced into the variant proteins, or the section 
    of the variants you want replaced with the constant protein. Each sequence of interest should be in their own 
    fasta file. All non-reference variant proteins will have the homologous sequence from an alignment replaces'''
    full_reference_fasta: str
    '''The path to full-length protein of the constant protein or the reference for the variants'''
    msa_file_name: str
    '''Path to the multiple sequence alignment to be created or a user-generated msa'''
    email_for_accession: str
    msa_fasta: str
    '''Name of the new fasta that will contain all variant protein sequences to be aligned'''
    constant_fasta_for_alphafold: str
    '''Path to the fasta to be predicted by Alphafold'''
    fasta_file_list_name: str
    '''filename to hold the complete list of fastas created in the fasta operation'''

    def __init__(self, fasta_dict):
        for key, value in fasta_dict.items():
            setattr(self, key, value)


class CrossSubmissionArguments:
    submission_toggles: dict
    '''A dictionary that holds optional operations within the submission operation, they are generally turned on by 
    placing True within the quotes of their value'''
    create_slurms: str
    sbatch_slurms: str
    stragglers_or_custom_or_all: str
    '''If stragglers is placed in the quotes it will check to see if an alphafold prediction is done for all created 
    chimeras and making a new the list of incomplete predictions or 'stragglers' to either create a file with their 
    fastas or put them into a slurm to be submitted again. If custom is placed, all chimeras within the 
    custom_list_to_run file will be submitted as slurm jobs. If all is selected all chimeras will be run 
    indiscriminately'''
    custom_list_to_run: str
    '''List that either be created by toggling on create_file_of_stragglers, or user-generated and should contain 
    line separated fastas of proteins you want predicted'''
    create_file_of_stragglers: str
    proteins_per_slurm: int
    template_slurm: str
    '''Path to the file that has the template settings for alphafold slurm jobs'''
    alphafold_shell_script: str
    '''Path to the shell script that runs alphafold, it must take a list of comma-separated fastas 
    and an output as its arguments. It's recommended you use the on provided in github'''
    slurm_output: str
    slurm_error: str

    def __init__(self, fasta_dict):
        for key, value in fasta_dict.items():
            setattr(self, key, value)


class CrossAnalysisArguments:
    analysis_toggles: dict
    make_plddts: str
    make_pdbs: str

    analysis_output_file: str
    '''Path of the file you want created that will contain all analysis created.'''
    column_names: dict
    '''A dictionary containing a multiple keys and their list value that each correspond with a type of data that is 
    output by the script, a user-selected title for the data column and whether that data is included in the final 
    output'''
    filename_stems: list
    '''A list ["True","Protein"] that contains quotes that when filled with True, include a data column in the final 
    output that contains the nicknames created from the naming convention previously established. The second index 
    controls the column name in the final output'''
    relative_stability: list
    overall_chimera_stability: list

    def __init__(self, fasta_dict):
        for key, value in fasta_dict.items():
            setattr(self, key, value)

# TODO add column names either class or method for gathering
class CrossChimeraContainer:
    argument_dict: dict
    operation_toggles: dict

    def __init__(self, shifted_json_file):
        with open(shifted_json_file, 'rb') as jfile:
            self.argument_dict = load(jfile)
        self.operation_toggles = self.argument_dict['operation_toggles']
        self.chimeras = ()

    def get_dict_args(self, dict_class, attr_name, arg_key):
        dict_args = self.argument_dict[arg_key]
        setattr(self, attr_name, dict_class(dict_args))

    def add_chimera(self, chimera):
        self.chimeras += (chimera,)
        
        
container_of_containers = ()
for index, jsn in enumerate(args.arg_jsons.split(',')):
    container=CrossChimeraContainer(jsn)
    container_of_containers += (container,)
    container.get_dict_args(CrossFastaArguments, 'fasta_args', 'fasta_arguments')
    if container.fasta_args.fasta_toggles['constant_or_variant']=='constant':
        constant_container=container
        constant_seq_of_interest = container.fasta_args.sequence_of_interest
        constant_fasta=container.fasta_args.full_reference_fasta
        constant_fasta_identifier = container.fasta_args.fasta_identifier

    elif container.fasta_args.fasta_toggles['constant_or_variant']=='variant':
        variant_container=container
        variant_seq_of_interest = container.fasta_args.sequence_of_interest
        variant_fasta=container.fasta_args.full_reference_fasta
        variant_fasta_identifier=container.fasta_args.fasta_identifier
fasta_toggles = variant_container.fasta_args.fasta_toggles


with open(constant_seq_of_interest, 'r') as fasta:
    constant_seq_of_interest = ''.join(x for x in fasta if x[0] != '>' if x != '').strip().replace('\n', '')
with open(variant_seq_of_interest, 'r') as fasta:
    variant_seq_of_interest = ''.join(x for x in fasta if x[0] != '>' if x != '').strip().replace('\n', '')

variant_container.get_dict_args(CrossFastaArguments, 'naming_args', 'naming_arguments')
placeholder = variant_container.naming_args.placeholder
msa=variant_container.fasta_args.msa_file_name
subunits = variant_container.fasta_args.number_of_subunits
alphafold_dir=variant_container.naming_args.alphafold_outputs_dir

if path.exists(variant_container.fasta_args.msa_file_name):
    with open(variant_container.fasta_args.msa_file_name, "r") as alignment:
        alignment = alignment.read().split('>')
        sequence_dictionary = {sequence.split('\n')[0]: ''.join(sequence.split('\n')[1:]) for sequence in alignment
                               if
                               len(sequence) != 0}
        for fasta_id, sequence in sequence_dictionary.items():
            chimera = chimeracls()
            variant_container.add_chimera(chimera)
            chimera.nickname=fasta_id
            chimera.native_seq=sequence


else:
    with open(variant_container.fasta_args.protein_list, 'r') as info_list:
        info_list = info_list.readlines()
        for line in info_list:
            chimera=chimeracls()
            variant_container.add_chimera(chimera)
            chimera.nickname=line.split()[1]
            chimera.accession = line.split()[0]
list_of_chis=variant_container.chimeras
num_of_chi=len(list_of_chis)
for chimera in list_of_chis:
    chimera.monomer_stem=variant_container.naming_args.monomer_naming_convention.replace(placeholder, chimera.nickname)
    chimera.chimera_stem=variant_container.naming_args.chimera_naming_convention.replace(placeholder, chimera.nickname)
    chimera.chi_pdb = path.join(f'{alphafold_dir}{chimera.chimera_stem}', 'ranked_0.pdb')
    if subunits>1:
        #make multimer = monomer when subunits is greater
        chimera.multimer_stem=variant_container.naming_args.multimer_naming_convention.replace(placeholder, chimera.nickname)
        chimera.multimer_fasta = variant_container.naming_args.fasta_directory + chimera.multimer_stem + variant_container.naming_args.fasta_extension
        chimera.native_pdb= path.join(f'{alphafold_dir}{chimera.multimer_stem}', 'ranked_0.pdb')
    else:
        chimera.native_pdb = path.join(f'{alphafold_dir}{chimera.monomer_stem}', 'ranked_0.pdb')
if variant_container.operation_toggles['run_fasta_operation']:

    for chimera in list_of_chis:
        chimera.monomer_fasta = variant_container.naming_args.fasta_directory + chimera.monomer_stem + variant_container.naming_args.fasta_extension
        chimera.chimera_fasta = variant_container.naming_args.fasta_directory + chimera.chimera_stem + variant_container.naming_args.fasta_extension
        # use sequences from alignment if user-generated
    if fasta_toggles['Create_an_alignment'] or not path.exists(variant_container.fasta_args.msa_file_name):
        email = variant_container.fasta_args.email_for_accession
        for chimera in list_of_chis:
            if subunits == 1:
                accession_to_fasta(chimera.monomer_fasta, chimera.accession, email, subunits)
            else:
                accession_to_fasta(chimera.monomer_fasta, chimera.accession, email, subunits,
                          chimera.multimer_fasta)
        multiple_sequence_alignment(tuple(chimera.monomer_fasta for chimera in list_of_chis), variant_container.fasta_args.msa_fasta, msa,
                                                         variant_fasta, variant_container.fasta_args.muscle_command_for_msa)
    for chimera in list_of_chis:
        variant_splice=alignment_finder(msa, variant_seq_of_interest, chimera.nickname,
                                                               variant_fasta_identifier)[0]
        chimera.chi_seq=chimera.native_seq.replace(variant_splice,constant_seq_of_interest)
        if subunits == 1:
            fasta_creation(chimera.monomer_fasta,[tuple((chimera.native_seq,subunits,chimera.monomer_stem))])
        else:
            fasta_creation(chimera.multimer_fasta, [tuple((chimera.native_seq, subunits, chimera.multimer_stem))])
        fasta_creation(chimera.chimera_fasta, [tuple((chimera.chi_seq,subunits,chimera.chimera_stem))])
    if fasta_toggles['Make_a_list_of_created_fasta_files']:
        with open(variant_container.fasta_args.fasta_list_file_name, 'w') as fasta_list_file:
            fasta_list_file.write("\n".join(chimera.chimera_fasta for chimera in list_of_chis)+'\n')
            if subunits == 1:
                fasta_list_file.write("\n".join(chimera.monomer_fasta for chimera in list_of_chis)+'\n')
            else:
                fasta_list_file.write("\n".join(chimera.multimer_fasta for chimera in list_of_chis)+'\n')
    
# if operation_toggles['run_analysis_operation']:
#     plddt_directory,plddt_extension = variant_container.naming_args.plddt_directory,variant_container.naming_args.plddt_extension
#     analysis_toggles = analysis_arguments.analysis_toggles
#     constant_plddt = analysis_arguments.constant_plddt
#     chimera_boundary = [(0, len(constant_seq_of_interest)),
#                         (len(constant_seq_of_interest), None)]
#     if prefix_protein == 'variant': chimera_boundary = chimera_boundary[::-1]
#     averaged_constant_plddt = analysis_arguments.averaged_constant
#     constant_boundary = \
#         AccessiontoAlignment.alignment_finder(constant_fasta, constant_seq_of_interest, constant_fasta_identifier,
#                                               constant_fasta_identifier)[
#             1]
#     if analysis_toggles.make_plddts:
#         Analysis.generate_alphafold_files(
#             analysis_arguments.alphafold_outputs_directory + analysis_arguments.constant_alphafold_folder_name + '/',
#             constant_plddt)
#         if subunits > 1:
#             Analysis.averaging_multimer_plddt(constant_plddt,
#                                               analysis_arguments.averaged_constant,
#                                               subunits)
#     for chimera in list_of_chis:
#         chimera.native_plddt=plddt_directory + chimera.monomer_stem + plddt_extension if subunits==1 else plddt_directory + chimera.multimer_stem + plddt_extension
#         chimera.chi_plddt=plddt_directory + chimera.chimera_stem + plddt_extension
#         native_folder=analysis_arguments.alphafold_outputs_directory + chimera.monomer_stem + '/' if subunits==1 else analysis_arguments.alphafold_outputs_directory + chimera.multimer_stem + '/'
#         chimera_folder=analysis_arguments.alphafold_outputs_directory + chimera.chimera_stem + '/'
#         if analysis_toggles.make_pdbs:
#             Analysis.generate_alphafold_files(native_folder,"NA",chimera.native_pdb)
#             Analysis.generate_alphafold_files(chimera_folder, "NA", chimera.chi_pdb)
#         if subunits>1:
#             avg_native_plddt = plddt_directory + variant_container.naming_args.averaged_multimer_naming_convention.replace(placeholder, chimera.nickname) + plddt_extension
#             avg_chimera_plddt = plddt_directory + variant_container.naming_args.averaged_chimera_naming_convention.replace(placeholder,
#                                                                                                  chimera.nickname) + plddt_extension
#         if analysis_toggles['make_plddts']:
#             Analysis.generate_alphafold_files(native_folder, chimera.native_plddt)
#             Analysis.generate_alphafold_files(chimera_folder, chimera.chi_plddt)
#             if subunits > 1:
#                 Analysis.averaging_multimer_plddt(chimera.native_plddt, avg_native_plddt, subunits )
#                 Analysis.averaging_multimer_plddt(chimera.chi_plddt, avg_chimera_plddt, subunits)
#         if analysis_toggles['make_emboss_files']:
#             emboss=variant_container.naming_args.emboss_names.replace(placeholder, chimera.nickname)
#             variant_splice_info =AccessiontoAlignment.alignment_finder(msa, variant_seq_of_interest, chimera.nickname,
#                                                                        variant_fasta_identifier, analysis_arguments.emboss_command,
#                                                                        emboss)[1]
#         else:
#             variant_splice_info = AccessiontoAlignment.alignment_finder(msa, variant_seq_of_interest, chimera.nickname,
#                                                                         variant_fasta_identifier)[1]
#         if subunits == 1:
#             constant_stability = Analysis.relative_stability(constant_plddt, constant_boundary,
#                                                              chimera.chi_plddt, chimera_boundary[0])
#             variant_stability = Analysis.relative_stability(chimera.native_plddt, variant_splice_info,
#                                                             chimera.chi_plddt, chimera_boundary[1])
#             chimera.rel_stability=(constant_stability[0] + variant_stability[0]) / (constant_stability[1] + variant_stability[1])
#         else:
#             constant_stability = Analysis.relative_stability(averaged_constant_plddt, constant_boundary,
#                                                              avg_chimera_plddt,
#                                                              chimera_boundary[0])
#             variant_stability = Analysis.relative_stability(avg_native_plddt, variant_splice_info,
#                                                             avg_chimera_plddt, chimera_boundary[1])
#             chimera.rel_stability=(constant_stability[0] + variant_stability[0]) / (constant_stability[1] + variant_stability[1])
#     column_names = analysis_arguments.column_names
#     data_dict= {'average_relative_stability': tuple(chimera.rel_stability for chimera in list_of_chis),'nickname': tuple(chimera.nickname for chimera in list_of_chis)}
#     if column_names['similarity'][0]:
#         emboss_files = tuple(variant_container.naming_args.emboss_names.replace(placeholder, chimera.nickname) for chimera in
#                              list_of_chis)
#         data_dict.similarity = tuple(map(Analysis.get_sequence_similarity, emboss_files))
#     if column_names.overall_native_stability[0]:
#         data_dict.overall_native_stability = tuple(Analysis.overall_confidence(chimera.native_plddt) for chimera in list_of_chis)
#     if column_names.overall_chimera_stability[0]:
#         data_dict.overall_chimera_stability = tuple(Analysis.overall_confidence(chimera.chi_plddt) for chimera in list_of_chis)
#     column_names = {key: value[1] for (key, value) in analysis_arguments.column_names.items() if value[0]}
#     data_array = empty((len(data_dict.average_relative_stability) + 1, len(column_names)), dtype=object)
#     for column_count, (corresponding_data, column_name) in enumerate(column_names.items()):
#         data_array[0, column_count], data_array[1:, column_count] = column_name, data_dict[corresponding_data]
#     savetxt(analysis_arguments.analysis_output_csv, data_array, fmt=','.join('%s' for x in column_names),
#             delimiter=",")
#     del data_dict
#     # TODO have a functiuon to proint pdbs and make pdbs for reference
#     # TODO make it to be able to simulate native pdb
# # TODO: Create a way for python to check periodically whether the setup slurms are done and then submit the prodcution
# # TODO make certain parts available from commandline
# if operation_toggles['run_gromacs_operation']=='True':
#     gromacs_data_dict={}
#     gromacs_arguments = argument_dict.gromacs_arguments
#     gromacs_toggles = gromacs_arguments.gromacs_toggles
#     # TODO if you hyand select make a file to resubmit with the choisces
#     if not path.exists(gromacs_arguments.pdbs_to_run) or gromacs_toggles['create new pdb list']:
#         with open(analysis_arguments.analysis_output_csv) as data:
#             pdbs_to_run=tuple(chimera.chi_pdb for chimera in list_of_chis)
#             for row in reader(data):
#                 print(row)
#             for index,pdb in enumerate(pdbs_to_run):
#                 print(f'{index}. {pdb}')
#             pdbs_to_run=tuple(pdbs_to_run[int(index)] for index in input("Submit comma separated list:\n Ex: 0,1,2,7,8\n Enter:").split(','))
#             with open(gromacs_arguments.pdbs_to_run,'w') as pdb_list:
#                 for pdb in pdbs_to_run:
#                     pdb_list.write(f'{pdb}\n')
#     else:
#         with open(gromacs_arguments.pdbs_to_run,'r') as pdbs_to_run:
#             pdbs_to_run=tuple(pdb.split()[0] for pdb in pdbs_to_run.readlines())
#     for pdb in pdbs_to_run:
#         system(f'cp {pdb} {variant_container.naming_args["gromacs_slurm_dir"]}')
#     gromacs_data_dict.setup_slurms = tuple(
#         variant_container.naming_args.gromacs_slurm_dir + Path(pdb).stem + variant_container.naming_args.setup_extension for pdb in
#         pdbs_to_run)
#     # TODO make folders for gromacs
#     if gromacs_toggles['create setup slurms']=='True':
#         for index,(pdb,slurm) in enumerate(zip(pdbs_to_run,gromacs_data_dict.setup_slurms)):
#             create_setup_slurm(pdb,gromacs_arguments.gmxbin,gromacs_arguments.pdb2gmx,gromacs_arguments.slurm_template,slurm,gromacs_arguments.slurm_output.replace(placeholder,Path(pdb).stem),gromacs_arguments.slurm_error.replace(placeholder,Path(pdb).stem))
#     if gromacs_toggles['sbatch setup slurms']:
#         for slurm in gromacs_data_dict.setup_slurms:
#             system(f"sbatch {slurm}")
#     gromacs_data_dict.mdrun_slurms = tuple(
#         variant_container.naming_args.gromacs_slurm_dir + Path(pdb).stem + variant_container.naming_args.production_extension for pdb
#         in pdbs_to_run)
#     if gromacs_toggles['create mdrun slurms']:
#         for index, (pdb, slurm) in enumerate(zip(pdbs_to_run, gromacs_data_dict.mdrun_slurms)):
#             create_prod_slurm(pdb,gromacs_arguments.gmxbin,gromacs_arguments.slurm_template,slurm,gromacs_arguments.slurm_output.replace(placeholder,Path(pdb).stem),gromacs_arguments.slurm_error.replace(placeholder,Path(pdb).stem))
#     if gromacs_toggles['sbatch mdrun slurms']:
#         for slurm in gromacs_data_dict.mdrun_slurms:
#             system(f"sbatch {slurm}")
