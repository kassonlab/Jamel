import copy
from json import load
from os import system, path
from pathlib import Path
from sys import exit
from time import perf_counter
from numpy import empty, savetxt
import Analysis
from setup import create_alphafold_slurm
from AccessiontoAlignment import alignment_finder, accession_to_fasta, multiple_sequence_alignment
from ChimeraGenerator import fasta_creation, chimeracls, update_json
from itertools import product
import argparse

# TODO add a commandline to change specific keys and replace recurring string through out all keys
# and be able to request operation from command line that overwrites whats in json
parser = argparse.ArgumentParser(
    description='Creating chimeric proteins from one or more json proteins corresponding to each unique protein in '
                'the structure')
subparser = parser.add_subparsers(dest='cmdline_toggles')
parser.add_argument('-u', '--updatejson', type=str, required=False,
                    help='updating json configs ex. default_json,old_json')
parser.add_argument('-i', '--jsoninput', dest='arg_jsons', required=False, type=str,
                    help='Json config file input, comma separated for each unique protein. Make sure initial json has '
                         'naming conventions')
operations = subparser.add_parser('operations')
operations.add_argument('-fa', '--fasta', required=False, action='store_true', default=False,
                        help='turns on fasta operation')
operations.add_argument('-s', '--submission', required=False, action='store_true', default=False,
                        help='turns on submission operation')
operations.add_argument('-a', '--analysis', required=False, action='store_true', default=False,
                        help='turns on analysis operation')
args = parser.parse_args()
operation_args = (args.fasta,args.submission,args.analysis)
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
    placing true within the quotes of their value'''
    Make_a_lst_of_created_fasta_files: bool
    '''Create a file with line separted list of all fastas created'''
    Create_an_alignment: bool
    constant_or_variant: bool
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
    placing true within the quotes of their value'''
    create_slurms: bool
    sbatch_slurms: bool
    stragglers_or_custom_or_all: bool
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
    make_plddts: bool
    make_pdbs: bool
    make_emboss_files:bool

    emboss_command: str
    analysis_output_csv: str
    '''Path of the file you want created that will contain all analysis created.'''
    column_names: dict
    '''A dictionary containing a multiple keys and their list value that each correspond with a type of data that is 
    output by the script, a user-selected title for the data column and whether that data is included in the final 
    output'''
    similarity: list
    overall_native_stability: list
    nickname: list
    '''A list [true,"Protein"] that contains quotes that when filled with true, include a data column in the final 
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
    run_fasta_operation: bool
    alphafold_submission: bool
    run_analysis_operation: bool
    run_gromacs_operation: bool

    def __init__(self, shifted_json_file):
        self.naming_args = None
        self.submission_args = None
        self.fasta_args = None
        self.naming_args = None
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
    container = CrossChimeraContainer(jsn)
    container_of_containers += (container,)
    container.get_dict_args(CrossFastaArguments, 'fasta_args', 'fasta_arguments')
    if container.fasta_args.fasta_toggles['constant_or_variant'] == 'constant':
        constant_container = container
        constant_seq_of_interest = container.fasta_args.sequence_of_interest
        constant_fasta = container.fasta_args.full_reference_fasta
        constant_fasta_identifier = container.fasta_args.fasta_identifier
        constant_submission = container.fasta_args.constant_fasta_for_alphafold
    elif container.fasta_args.fasta_toggles['constant_or_variant'] == 'variant':
        variant_container = container
        variant_seq_of_interest = container.fasta_args.sequence_of_interest
        variant_fasta = container.fasta_args.full_reference_fasta
        variant_fasta_identifier = container.fasta_args.fasta_identifier
fasta_toggles = variant_container.fasta_args.fasta_toggles
if any(operation_args):
    for key in variant_container.operation_toggles:
        variant_container.operation_toggles[key] = False

with open(constant_seq_of_interest, 'r') as fasta:
    constant_seq_of_interest = ''.join(x for x in fasta if x[0] != '>' if x != '').strip().replace('\n', '')
with open(variant_seq_of_interest, 'r') as fasta:
    variant_seq_of_interest = ''.join(x for x in fasta if x[0] != '>' if x != '').strip().replace('\n', '')

variant_container.get_dict_args(CrossFastaArguments, 'naming_args', 'naming_arguments')
placeholder = variant_container.naming_args.placeholder
msa = variant_container.fasta_args.msa_file_name
subunits = variant_container.fasta_args.number_of_subunits
alphafold_dir = variant_container.naming_args.alphafold_outputs_dir

if path.exists(variant_container.fasta_args.msa_file_name):
    with open(variant_container.fasta_args.msa_file_name, "r") as alignment:
        alignment = alignment.read().split('>')
        sequence_dictionary = {sequence.split('\n')[0]: ''.join(sequence.split('\n')[1:]) for sequence in alignment
                               if
                               len(sequence) != 0}
        for fasta_id, sequence in sequence_dictionary.items():
            chimera = chimeracls()
            variant_container.add_chimera(chimera)
            chimera.nickname = fasta_id
            chimera.native_seq = sequence.replace('-', '')
else:
    with open(variant_container.fasta_args.protein_list, 'r') as info_list:
        info_list = info_list.readlines()
        for line in info_list:
            chimera = chimeracls()
            variant_container.add_chimera(chimera)
            chimera.nickname = line.split()[1]
            chimera.accession = line.split()[0]

list_of_chis = variant_container.chimeras
num_of_chi = len(list_of_chis)
for chimera in list_of_chis:
    chimera.monomer_stem = variant_container.naming_args.monomer_naming_convention.replace(placeholder,
                                                                                           chimera.nickname)
    chimera.chimera_stem = variant_container.naming_args.chimera_naming_convention.replace(placeholder,
                                                                                           chimera.nickname)
    chimera.chi_pdb = path.join(f'{alphafold_dir}{chimera.chimera_stem}', 'ranked_0.pdb')
    chimera.monomer_fasta = variant_container.naming_args.fasta_directory + chimera.monomer_stem + variant_container.naming_args.fasta_extension
    chimera.chimera_fasta = variant_container.naming_args.fasta_directory + chimera.chimera_stem + variant_container.naming_args.fasta_extension
    # make multimer = monomer when subunits is greater
    chimera.multimer_stem = variant_container.naming_args.multimer_naming_convention.replace(placeholder,
                                                                                             chimera.nickname)
    chimera.multimer_fasta = variant_container.naming_args.fasta_directory + chimera.multimer_stem + variant_container.naming_args.fasta_extension

    if subunits == 1:
        chimera.multimer_stem = variant_container.naming_args.monomer_naming_convention.replace(placeholder,
                                                                                                 chimera.nickname)
        chimera.multimer_fasta = variant_container.naming_args.fasta_directory + chimera.monomer_stem + variant_container.naming_args.fasta_extension

    chimera.native_pdb = path.join(f'{alphafold_dir}{chimera.multimer_stem}', 'ranked_0.pdb')
    chimera.chi_pdb = path.join(f'{alphafold_dir}{chimera.chimera_stem}', 'ranked_0.pdb')

if variant_container.operation_toggles['run_fasta_operation'] or args.fasta:
    if fasta_toggles['Create_an_alignment'] or not path.exists(variant_container.fasta_args.msa_file_name):
        email = variant_container.fasta_args.email_for_accession
        for chimera in list_of_chis:
            accession_to_fasta(chimera.monomer_fasta, chimera.accession, email, subunits,
                               chimera.multimer_fasta)
        multiple_sequence_alignment(tuple(chimera.monomer_fasta for chimera in list_of_chis),
                                    variant_container.fasta_args.msa_fasta, msa,
                                    variant_fasta, variant_container.fasta_args.muscle_command_for_msa)
    for chimera in list_of_chis:
        variant_splice = alignment_finder(msa, variant_seq_of_interest, chimera.nickname,
                                          variant_fasta_identifier)[0]
        chimera.chi_seq = chimera.native_seq.replace(variant_splice, constant_seq_of_interest)
        fasta_creation(chimera.multimer_fasta, [tuple((chimera.native_seq, subunits, chimera.multimer_stem))])
        fasta_creation(chimera.chimera_fasta, [tuple((chimera.chi_seq, subunits, chimera.chimera_stem))])
    if fasta_toggles['Make_a_list_of_created_fasta_files']:
        with open(variant_container.fasta_args.fasta_list_file_name, 'w') as fasta_list_file:
            fasta_list_file.write(constant_submission + '\n')
            fasta_list_file.write("\n".join(chimera.chimera_fasta for chimera in list_of_chis) + '\n')
            fasta_list_file.write("\n".join(chimera.multimer_fasta for chimera in list_of_chis) + '\n')

if variant_container.operation_toggles['alphafold_submission'] or args.submission:
    variant_container.get_dict_args(CrossSubmissionArguments, 'submission_args', 'alphafold_submission_args')
    submission_toggles = variant_container.submission_args.submission_toggles
    output_directory = variant_container.naming_args.alphafold_outputs_dir
    fastas = [chimera.chimera_fasta for chimera in list_of_chis] + [constant_submission]
    fastas=fastas + [chimera.multimer_fasta for chimera in list_of_chis]
    fasta_to_run = ()
    proteins_per_slurm = variant_container.submission_args.proteins_per_slurm
    template_slurm = variant_container.submission_args.template_slurm
    alphafold_shell_script = variant_container.submission_args.alphafold_shell_script
    naming_convention = variant_container.submission_args.slurm_naming
    # Loops through all fastas created and checks if they are complete by looking for their ranking_debug file
    if submission_toggles['stragglers_or_custom_or_all'] == 'stragglers':
        for fasta in fastas:
            if not path.exists(output_directory + Path(fasta).stem + '/ranking_debug.json'):
                fasta_to_run += (fasta,)
        # Puts all fastas in a line separated file specified by custom_list_to_run
        if submission_toggles['create file of stragglers']:
            with open(variant_container.submission_args.custom_list_to_run, 'w') as run_list:
                run_list.write('\n'.join(fasta for fasta in fasta_to_run))
        # if all of them are complete and fasta_to_run is empty then all slurm actions are toggled off
        if not fasta_to_run:
            submission_toggles['create_slurms'] = submission_toggles['sbatch slurms'] = ''
    elif submission_toggles['stragglers_or_custom_or_all'] == 'custom':
        with open(variant_container.submission_args.custom_list_to_run, 'r') as run_list:
            run_list = run_list.readlines()
            fasta_to_run = [x.split()[0] for x in run_list]
    elif submission_toggles['stragglers_or_custom_or_all'] == 'all':
        fasta_to_run = fastas
    if submission_toggles['create_slurms']:
        for slurm_index, file_index in enumerate(range(0, len(fasta_to_run), proteins_per_slurm)):
            current_slurm = naming_convention.replace(placeholder, str(slurm_index))
            create_alphafold_slurm(fasta_to_run[file_index:file_index + proteins_per_slurm], current_slurm,
                                   template_slurm,
                                   variant_container.submission_args.slurm_output.replace(placeholder,
                                                                                          str(slurm_index)),
                                   variant_container.submission_args.slurm_error.replace(placeholder,
                                                                                         str(slurm_index)),
                                   alphafold_shell_script, output_directory)
    if submission_toggles['sbatch slurms']:
        for slurm_index, file_index in enumerate(range(0, len(fasta_to_run), proteins_per_slurm)):
            current_slurm = naming_convention.replace(placeholder, str(slurm_index))
            system(f'sbatch {current_slurm}')

if variant_container.operation_toggles['run_analysis_operation'] or args.analysis:
    alphafold_directory = variant_container.naming_args.alphafold_outputs_dir
    plddt_direc, plddt_ext, pdb_direc, pdb_ext = variant_container.naming_args.plddt_directory, \
        variant_container.naming_args.plddt_extension, variant_container.naming_args.pdb_directory, \
        variant_container.naming_args.pdb_extension
    reference_stem = Path(constant_submission).stem
    variant_container.get_dict_args(CrossAnalysisArguments, 'analysis_args', 'analysis_arguments')
    analysis_toggles = variant_container.analysis_args.analysis_toggles
    reference_pdb = path.join(f'{alphafold_directory}{reference_stem}', 'ranked_0.pdb')
    reference_plddt = Analysis.get_plddt_tuple_from_pdb(reference_pdb)

    for chimera in list_of_chis:
        chimera.native_plddt= Analysis.get_plddt_tuple_from_pdb(chimera.native_pdb)[chimera.native_seq]
        chimera.chi_plddt= Analysis.get_plddt_tuple_from_pdb(chimera.native_pdb)[chimera.chi_seq]

    if analysis_toggles['make_plddts']:
        for chimera in list_of_chis:
            Analysis.get_plddt_file_from_pdb(chimera.native_pdb,
                                             path.join(plddt_direc, chimera.multimer_stem + plddt_ext))
            Analysis.get_plddt_file_from_pdb(chimera.chi_pdb,
                                             path.join(plddt_direc, chimera.chimera_stem + plddt_ext))
        Analysis.get_plddt_file_from_pdb(reference_pdb, path.join(plddt_direc, reference_stem + plddt_ext))

    if analysis_toggles["make_pdbs"]:
        for chimera in list_of_chis:
            Analysis.generate_alphafold_files(f'{alphafold_directory}{chimera.multimer_stem}',
                                              new_pdb=path.join(pdb_direc, chimera.multimer_stem + pdb_ext))
            Analysis.generate_alphafold_files(f'{alphafold_directory}{chimera.chimera_stem}',
                                              new_pdb=path.join(pdb_direc, chimera.chimera_stem + pdb_ext))
        Analysis.generate_alphafold_files(
            f'{alphafold_directory}{reference_stem}',
            new_pdb=path.join(pdb_direc,
                              reference_stem + pdb_ext))

    for chimera in list_of_chis:
        #relative stability for constant versus variant
        chimera.rel_stability = Analysis.average_relative_stability_full_chimera(chimera.native_plddt,
                                                                                 chimera.partner_boundaries,
                                                                                 chimera.chi_plddt,
                                                                                 reference_plddt,
                                                                                 constant_seq_of_interest,
                                                                                 constant_submission,
                                                                                 constant_fasta_identifier)

    data_columns = {}
    # Checks which data columns are wanted by the user by looking for True in the first index of each array in
    # column_names from the analysis_args, Each container can have its own column preferences and every container
    # will have its own columns of data per data column requested
    for container in container_of_containers:
        container.get_dict_args(CrossAnalysisArguments, 'analysis_args', 'analysis_arguments')
        column_choices = container.analysis_args.column_names
        # TODO add this functionality as a mehtod
        if column_choices['filename_stems'][0]:
            data_columns[column_choices['filename_stems'][1]] = tuple(
                chimera.file_stem for chimera in container.chimeras)
        if column_choices['relative_stability'][0]:
            data_columns[column_choices['relative_stability'][1]] = tuple(
                chimera.rel_stability for chimera in container.chimeras)
        if column_choices['overall_chimera_stability'][0]:
            data_columns[column_choices['overall_chimera_stability'][1]] = tuple(
                Analysis.overall_confidence(chimera.plddt) for chimera in container.chimeras)

    data_array = empty(((num_of_chis + 1), len(data_columns)), dtype=object)
    for column_count, (column_name, data) in enumerate(data_columns.items()):
        data_array[0, column_count] = column_name
        data_array[1:, column_count] = data
    savetxt(variant_container.analysis_args.analysis_output_file, data_array, fmt=','.join('%s' for x in data_columns))

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
