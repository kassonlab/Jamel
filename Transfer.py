from os import listdir,path,system,chdir
import argparse
from shutil import rmtree
from Analysis import skeletonize_alphafold_folder

parser = argparse.ArgumentParser(description='')
parser.add_argument('-a', '--start', type=str, required=False, help='')
parser.add_argument('-s', '--to', type=str, required=True, help='')
parser.add_argument('-z', '--zip', type=str, required=False, help='')
parser.add_argument('--delete', type=bool,required=False, default=False,help='')
args = parser.parse_args()
folder_of_folders=[folder for folder in listdir(args.start)]
if args.start and args.to:
    for folder in folder_of_folders:
        if skeletonize_alphafold_folder(path.join(args.start, folder),path.join(args.to, folder)) and path.exists(path.join(args.start, folder)) and args.delete:
            rmtree(path.join(args.start, folder))
if args.zip:
    chdir(args.to+'../')
    system(f'zip -m -r {args.zip}.zip {args.zip}')

