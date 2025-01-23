import copy
from json import load
from os import system, path
from pathlib import Path
from sys import exit
from time import perf_counter
from numpy import empty, savetxt
import Analysis
from setup import create_alphafold_slurm
from AccessiontoAlignment import alignment_finder
import ChimeraClasses
from ChimeraGenerator import update_json
from itertools import product
import argparse
# TODO allow change of config options from commandline
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
if __name__ == '__main__':
    TOTAL_ARGS = ChimeraClasses.ShiftedChimeraArgs(args.arg_jsons)
    if hasattr(args, 'fasta') or hasattr(args, 'submission') or hasattr(args, 'analysis'):
        operation_args = (args.fasta, args.submission, args.analysis)
        if any(operation_args) and operation_args:
            for key in TOTAL_ARGS.operation_toggles:
                TOTAL_ARGS.operation_toggles[key] = False
    if TOTAL_ARGS.operation_toggles['run_fasta_operation']:
        TOTAL_ARGS.fasta_operations()
