from json import load
from os import path
from pathlib import Path
from sys import exit

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import AccessiontoAlignment, Analysis, ChimeraGenerator
from setup import alphafold_submission_for_chimera_container
import argparse

# TODO add autocomplete for changing keys?? readline
# TODO be able to designate multiple seq_of_interest
# TODO labeled_schema_aln.csv accession route
# TODO labeled_schema_aln.csv this between Homomer proteins
# TODO be able tos swap into the constant
#TODO turn parser into function??
parser = argparse.ArgumentParser(
    description='Creating chimeric proteins where a region from a constant protein is spliced into a list of '
                'homologus regions between a protein family or vice versa')
parser.add_argument('-u', '--updatejson', type=str, required=False,
                    help='updating json configs. ex: default_json,old_json')
parser.add_argument('-i', '--jsoninput', dest='arg_jsons', required=False, type=str,
                    help='Comma seperated json config inputs. ex: constant_json,variant_json')
parser.add_argument('-ch', '--change', required=False, type=str,
                    help='input json with options you wish to change. ex: json')
parser.add_argument('-nv', '--new_values', required=False, type=str,
                    help='input comma seperated option:new_value pairs, seperate each part of the pair with a colon. '
                         'ex: option:new_value,option:new_value')
parser.add_argument('--overwrite', action='store_true', required=False, default=False,
                    help='This flag will overwrite the input json files for change,findnreplace, and updatejson, '
                         'the default without the flag will create a new json named new+input_file_stem.')
parser.add_argument('-pr', '--prints', required=False, type=str,
                    help='print all changeable options from inputted json or print current value for inputted option '
                         'in inputted json. ex: json or json,alphafold_outputs_dir')
parser.add_argument('-fr', '--findnreplace', required=False, type=str,
                    help='Warning:It is recommended you do not use overwrite flag with this option. Replaces a '
                         'recurring between a find:replace pair input within the json input in change flag. ex: '
                         'string_to_find:replacing_string')

subparser = parser.add_subparsers(dest='cmdline_toggles')
operations = subparser.add_parser('operations')
operations.add_argument('-fa', '--fasta', required=False, action='store_true', default=False,
                        help='turns on fasta operation. Flag must be called at the end of the command and '
                             'preceded by \'operations\' ex: python file.py -i json1,json2 operations -fa . Can be '
                             'combined with -s')
operations.add_argument('-s', '--submission', required=False, action='store_true', default=False,
                        help='turns on submission operation. This flag must be called at the end of the command and '
                             'preceded by \'operations\' ex: python file.py -i json1,json2 operations -s . Can be '
                             'combined with -fa or -a')
operations.add_argument('-a', '--analysis', required=False, action='store_true', default=False,
                        help='turns on analysis operation. This flag must be called at the end of the command and '
                             'preceded by \'operations\' ex: python file.py -i json1,json2 operations -a . Can be '
                             'combined with -s')

args = parser.parse_args()
# All conditionals to check for flags related to manipulation or clarification of json config files: (u,ch,pr,nv,overwite,fr)
# If any of them are called the program will cease after completion or error

if args.updatejson:
    args.updatejson = args.updatejson.split(',')
    ChimeraGenerator.update_json(args.updatejson[0], args.updatejson[1], args.overwrite)
    exit()
if args.prints:
    print_args = args.prints.split(',')
    with open(print_args[0], 'rb') as jfile:
        selected_json = load(jfile)
    if len(print_args) == 1:
        ChimeraGenerator.print_keys(selected_json)
    else:
        ChimeraGenerator.print_keys(selected_json, print_args[1])
    exit()
if args.new_values and args.change:
    new_key_values = [pair.split(':') for pair in args.new_values.split(',')]
    new_key_values = {key: value for key, value in new_key_values}
    for key, value in new_key_values.items():
        ChimeraGenerator.change_json_value(args.change, key, value, args.overwrite)
    exit()
if args.findnreplace and args.change:
    pair = tuple(args.findnreplace.split(':'))
    ChimeraGenerator.change_json_value(args.change, overwrite=args.overwrite, find_n_replace_tuple=pair)
    exit()


# Allows for commandline operation flags to override selections made in the json config


class HomomerNamingArguments:
    file_stem_placeholder: str
    WT_fasta_convention: str
    chimera_fasta_convention: str
    output_directory: str

    # All non-nested keys from the json are set as attributes to be called later
    def __init__(self, naming_dict):
        for key, value in naming_dict.items():
            setattr(self, key, value)


class HomomerFastaArguments:
    Make_a_lst_of_created_fasta_files: bool
    '''Create a file with line separted list of all fastas created'''
    Create_an_alignment: bool
    number_of_subunits: int
    sequence_of_interest: str
    '''This is the path to the fasta containing the sequence of either the constant protein you want spliced into the 
    variant proteins, or the section of the variants you want replaced with the constant protein. Each sequence of 
    interest should be in their own fasta file. All non-reference variant proteins will have the Homomer sequence 
    from an alignment replaced'''
    msa_file_name: str
    '''Path to the multiple sequence alignment to be created or a user-generated msa'''
    email_for_accession: str
    fasta_file_list_name: str
    '''filename to hold the complete list of fastas created in the fasta operation'''

    def __init__(self, fasta_dict):
        for key, value in fasta_dict.items():
            setattr(self, key, value)


class HomomerSubmissionArguments:
    submission_toggles: dict
    '''A dictionary that holds optional operations within the submission operation, they are generally turned on by 
    placing true within the quotes of their value'''
    create_slurms: bool
    sbatch_slurms: bool
    stragglers_or_custom_or_all: bool = 'stragglers'
    '''If stragglers is placed in the quotes it will check to see if an alphafold prediction is done for all created 
    chimeras and making a new the list of incomplete predictions or 'stragglers' to either create a file with their 
    fastas or put them into a slurm to be submitted again. If custom is placed, all chimeras within the 
    custom_list_to_run file will be submitted as slurm jobs. If all is selected all chimeras will be run 
    indiscriminately'''
    create_file_of_stragglers: bool
    custom_list_to_run: str
    '''Path to file that for list that either be created by toggling on create_file_of_stragglers, or user-generated that should contain 
    line separated fastas of proteins you want predicted'''
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


class HomomerAnalysisArguments:
    analysis_toggles: dict
    make_plddts: bool
    make_pdbs: bool
    make_emboss_files: bool

    emboss_command: str
    analysis_output_csv: str
    '''Path of the file you want created that will contain all analysis created.'''
    column_names: dict
    '''A dictionary containing a multiple keys and their list value that each correspond with a type of data that is 
    output by the script, a user-selected title for the data column and whether that data is included in the final 
    output'''
    file_stem: list
    '''An array [true,"Protein"] that contains a boolean that will control whether the data column containing 
    filestems/nicknames created from the naming convention previously established will be included in the data csv.  
    The second index controls the column name in the csv'''
    overall_native_stability: list
    '''See file_stem docstring'''
    relative_stability: list
    '''See file_stem docstring'''
    overall_chimera_stability: list
    '''See file_stem docstring'''

    def __init__(self, fasta_dict):
        for key, value in fasta_dict.items():
            setattr(self, key, value)


class HomomerChimeraArgs:
    argument_dict: dict
    '''Dictionary containing all arguments for running the chimera script'''
    operation_toggles: dict
    '''Dictionary containing boolean toggles for all major operation:fasta creation, alphafold submission, 
    analysis of alphafold results, and gromacs simulation submission.'''
    run_fasta_operation: bool
    alphafold_submission: bool
    run_analysis_operation: bool
    run_gromacs_operation: bool

    def __init__(self, arg_json):
        self.naming_args = None
        self.submission_args = None
        self.fasta_args = None
        self.analysis_args = None
        with open(arg_json, 'rb') as jfile:
            self.argument_dict = load(jfile)
        self.operation_toggles = self.argument_dict['operation_toggles']
        self.chimeras = ()
        self.get_dict_args(HomomerFastaArguments, 'fasta_args', 'fasta_arguments')
        self.get_dict_args(HomomerNamingArguments, 'naming_args', 'naming_arguments')
        self.get_dict_args(HomomerSubmissionArguments, 'submission_args', 'submission_arguments')

    def get_dict_args(self, dict_class, attr_name, arg_key):
        dict_args = self.argument_dict[arg_key]
        setattr(self, attr_name, dict_class(dict_args))


def create_file_path(directory, file_stem, extension=''):
    return Path(directory).joinpath(file_stem, extension)

def msa_to_fastas(msa_file,output,subunits):
    for label, sequence in AccessiontoAlignment.create_dictionary_from_alignment(msa_file).items():
        fasta_file = create_file_path(output, label+'.fa')
        AccessiontoAlignment.fasta_creation(fasta_file,
                                            [SeqRecord(Seq(AccessiontoAlignment.no_gap_sequence_from_alignment(sequence)), id=label, description="")
                                             for x in range(subunits)])

if __name__ == '__main__':
    TOTAL_ARGS = HomomerChimeraArgs(args.arg_jsons)
    fasta_toggles = TOTAL_ARGS.fasta_args.fasta_toggles

    if hasattr(args, 'fasta') or hasattr(args, 'submission') or hasattr(args, 'analysis'):
        operation_args = (args.fasta, args.submission, args.analysis)
        if any(operation_args) and operation_args:
            for key in TOTAL_ARGS.operation_toggles:
                TOTAL_ARGS.operation_toggles[key] = False

    PLACEHOLDER = TOTAL_ARGS.naming_args.filestem_placeholder
    msa = TOTAL_ARGS.fasta_args.msa_file_name
    subunits = TOTAL_ARGS.fasta_args.number_of_subunits
    TOTAL_ARGS.naming_args.output_directory = Path(TOTAL_ARGS.naming_args.output_directory)


    if TOTAL_ARGS.operation_toggles['run_fasta_operation']:
        email = TOTAL_ARGS.fasta_args.email_for_accession
        # checking if all wild-type and non-reference fastas exist or that there creation was toggled


            # fasta_file = create_file_path(TOTAL_ARGS.naming_args, TOTAL_ARGS.naming_args.chimera_fasta_convention,
            #                               label,
            #                               '.fa')
            # AccessiontoAlignment.fasta_creation(fasta_file,
            #                                     [SeqRecord(Seq(sequence.replace('-', '')), id=label, description="")
            #                                      for x in range(TOTAL_ARGS.fasta_args.number_of_subunits)])

    # TODO easy msa creation

    if TOTAL_ARGS.operation_toggles['alphafold_submission']:
        TOTAL_ARGS.get_dict_args(HomomerSubmissionArguments, 'submission_args', 'alphafold_submission_args')
        alphafold_submission_for_chimera_container(TOTAL_ARGS, fastas)

    if TOTAL_ARGS.operation_toggles['run_analysis_operation']:
        seq_of_interest = AccessiontoAlignment.extract_seq_from_fasta(TOTAL_ARGS.fasta_args.sequence_of_interest)

        for chimera in list_of_chis:
            homologous_splice = AccessiontoAlignment.alignment_finder(msa, seq_of_interest, chimera.file_stem,
                                                                      TOTAL_ARGS.fasta_args.reference_identifier)[0]
            chimera.chi_seq = ref_sequence.replace(seq_of_interest, homologous_splice)
        TOTAL_ARGS.get_dict_args(HomomerAnalysisArguments, 'analysis_args', 'analysis_arguments')
        analysis_toggles = TOTAL_ARGS.analysis_args.analysis_toggles
        ref_file_stem = Path(TOTAL_ARGS.fasta_args.fasta_for_alphafold).stem
        ref_pdb = path.join(str(alphafold_dir),'3mer6VSB', 'ranked_0.pdb')
        ref_plddt = {ref_sequence: Analysis.get_plddt_dict_from_pdb(ref_pdb)[ref_sequence]}
        if analysis_toggles['make_emboss_files']:
            for chimera in list_of_chis:
                AccessiontoAlignment.alignment_finder(msa, seq_of_interest, chimera.file_stem,
                                                      TOTAL_ARGS.fasta_args.reference_identifier,
                                                      TOTAL_ARGS.analysis_args.emboss_command,
                                                      TOTAL_ARGS.naming_args.emboss_naming.replace(name_placeholder,
                                                                                                   chimera.file_stem))
        if analysis_toggles["make_pdbs"] or any(not path.exists(chimera.native_pdb) for chimera in list_of_chis) or any(
                not path.exists(chimera.chi_pdb) for chimera in list_of_chis):
            for chimera in list_of_chis:
                Analysis.generate_alphafold_files(f'{alphafold_dir}{chimera.multimer_stem}',
                                                  new_pdb=TOTAL_ARGS.naming_args.pdb_naming.replace(file_placeholder,
                                                                                                    chimera.multimer_stem))
                Analysis.generate_alphafold_files(f'{alphafold_dir}{chimera.chimera_stem}',
                                                  new_pdb=TOTAL_ARGS.naming_args.pdb_naming.replace(file_placeholder,
                                                                                                    chimera.chimera_stem))
            Analysis.generate_alphafold_files(
                f'{alphafold_dir}{ref_file_stem}',
                new_pdb=TOTAL_ARGS.naming_args.pdb_naming.replace(file_placeholder, ref_file_stem))
        for chimera in list_of_chis:
            chimera.native_plddt = {
                chimera.native_seq: Analysis.get_plddt_dict_from_pdb(chimera.native_pdb)[chimera.native_seq]}
            chimera.chi_plddt = {chimera.chi_seq: Analysis.get_plddt_dict_from_pdb(chimera.chi_pdb)[chimera.chi_seq]}
            chimera.overall_native_stability = Analysis.overall_confidence(chimera.native_plddt[chimera.native_seq])
            chimera.overall_chimera_stability = Analysis.overall_confidence(chimera.chi_plddt[chimera.chi_seq])
            chimera.sequence_similarity = Analysis.get_sequence_similarity(
                TOTAL_ARGS.naming_args.emboss_naming.replace(name_placeholder, chimera.file_stem))

        if analysis_toggles['make_plddts']:
            for chimera in list_of_chis:
                Analysis.create_plddt_file_from_pdb(chimera.native_pdb,
                                                    TOTAL_ARGS.naming_args.plddt_naming.replace(file_placeholder,
                                                                                             chimera.multimer_stem))
                Analysis.create_plddt_file_from_pdb(chimera.chi_pdb,
                                                    TOTAL_ARGS.naming_args.plddt_naming.replace(file_placeholder,
                                                                                             chimera.chimera_stem))
            Analysis.create_plddt_file_from_pdb(ref_pdb,
                                                TOTAL_ARGS.naming_args.plddt_naming.replace(file_placeholder,
                                                                                         ref_file_stem))
        for chimera in list_of_chis:
            homologous_splice = AccessiontoAlignment.alignment_finder(msa, seq_of_interest, chimera.file_stem,
                                                                      TOTAL_ARGS.fasta_args.reference_identifier)[0]
            chimera.rel_stability = Analysis.revamped_rs(chimera.native_plddt, chimera.chi_plddt, ref_plddt,
                                                         homologous_splice)

        data_columns = Analysis.determine_columns_from_container(TOTAL_ARGS)
        Analysis.convert_data_dict_to_csv(data_columns, TOTAL_ARGS)

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
    #         system(f'cp {pdb} {container.naming_args["gromacs_slurm_dir"]}')
    #     gromacs_data_dict.setup_slurms = tuple(
    #         container.naming_args.gromacs_slurm_dir + Path(pdb).stem + container.naming_args.setup_extension for pdb in
    #         pdbs_to_run)
    #     # TODO make folders for gromacs
    #     if gromacs_toggles['create setup slurms']=='True':
    #         for index,(pdb,slurm) in enumerate(zip(pdbs_to_run,gromacs_data_dict.setup_slurms)):
    #             create_setup_slurm(pdb,gromacs_arguments.gmxbin,gromacs_arguments.pdb2gmx,gromacs_arguments.slurm_template,slurm,gromacs_arguments.slurm_output.replace(placeholder,Path(pdb).stem),gromacs_arguments.slurm_error.replace(placeholder,Path(pdb).stem))
    #     if gromacs_toggles['sbatch setup slurms']:
    #         for slurm in gromacs_data_dict.setup_slurms:
    #             system(f"sbatch {slurm}")
    #     gromacs_data_dict.mdrun_slurms = tuple(
    #         container.naming_args.gromacs_slurm_dir + Path(pdb).stem + container.naming_args.production_extension for pdb
    #         in pdbs_to_run)
    #     if gromacs_toggles['create mdrun slurms']:
    #         for index, (pdb, slurm) in enumerate(zip(pdbs_to_run, gromacs_data_dict.mdrun_slurms)):
    #             create_prod_slurm(pdb,gromacs_arguments.gmxbin,gromacs_arguments.slurm_template,slurm,gromacs_arguments.slurm_output.replace(placeholder,Path(pdb).stem),gromacs_arguments.slurm_error.replace(placeholder,Path(pdb).stem))
    #     if gromacs_toggles['sbatch mdrun slurms']:
    #         for slurm in gromacs_data_dict.mdrun_slurms:
    #             system(f"sbatch {slurm}")
