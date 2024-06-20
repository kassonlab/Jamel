from os import listdir,path,system,chdir
import argparse
from shutil import rmtree
from Chimeragenesis.Analysis import skeletonize_alphafold_folder

parser = argparse.ArgumentParser(description='')
parser.add_argument('-s', '--start', type=str, required=True, help='Starting directory holding alphafold data')
parser.add_argument('-t', '--to', type=str, required=False, help='Storage directory')
parser.add_argument('-z', '--zip', type=str, required=False, help='Name of the zip file (including directory if '
                                                                  'necessary)')
parser.add_argument('--delete', action='store_true',required=False, default=False,help='Do you want to delete the folders once tarnsfer is complete?')
args = parser.parse_args()
folders_of_start=[folder for folder in listdir(args.start)]
if args.start and args.to:
    for folder in folders_of_start:
        if skeletonize_alphafold_folder(path.join(args.start, folder),path.join(args.to, folder)) and path.exists(path.join(args.start, folder)) and args.delete:
            rmtree(path.join(args.start, folder))
if args.zip and args.to and args.start:
    chdir(args.to+'../')
    system(f'zip -m -r {args.zip+".zip"} {args.to}')

