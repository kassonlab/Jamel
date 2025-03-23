from sys import exit
import ChimeraClasses
from ChimeraGenerator import update_json
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
    update_json(*args.updatejson.split(','))
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
    if TOTAL_ARGS.operation_toggles['alphafold_submission']:
        TOTAL_ARGS.submission_operations()
    if TOTAL_ARGS.operation_toggles['run_analysis_operation']:
        TOTAL_ARGS.analysis_operations()

