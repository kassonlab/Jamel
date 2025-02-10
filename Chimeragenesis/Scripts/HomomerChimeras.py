import os
from os import path
from pathlib import Path
from sys import exit
from json import load

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from ChimeraClasses import HomomerChimeraArgs,HomomerAnalysisArguments,HomomerFastaArguments,HomomerNamingArguments,HomomerSubmissionArguments
import AccessiontoAlignment, Analysis, ChimeraGenerator
import ESM
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
#python C:\Users\jamel\PycharmProjects\Jamel\Chimeragenesis\HomomerChimeras.py -i C:\Users\jamel\PycharmProjects\Jamel\Chimeragenesis\homologous.json

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


if __name__ == '__main__':
    TOTAL_ARGS = HomomerChimeraArgs(args.arg_jsons)
    if hasattr(args, 'fasta') or hasattr(args, 'submission') or hasattr(args, 'analysis'):
        operation_args = (args.fasta, args.submission, args.analysis)
        if any(operation_args) and operation_args:
            for key in TOTAL_ARGS.operation_toggles:
                TOTAL_ARGS.operation_toggles[key] = False
    if TOTAL_ARGS.operation_toggles['run_fasta_operation']:
        TOTAL_ARGS.fasta_operations()
    # TODO easy msa creation

    if TOTAL_ARGS.operation_toggles['alphafold_submission']:
        alphafold_submission_for_chimera_container(TOTAL_ARGS, TOTAL_ARGS.all_fastas)

    if TOTAL_ARGS.operation_toggles['run_analysis_operation']:
        TOTAL_ARGS.analysis_operations()

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
