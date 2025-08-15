from os import system
from pathlib import Path
import argparse
from shutil import rmtree, copy, copytree
import ChimeraClasses
from Analysis import skeletonize_alphafold_folder, skeletionize_gmx_folder

parser = argparse.ArgumentParser(description='')
parser.add_argument('-s', '--start', type=str, required=True, help='Starting directory holding data')
parser.add_argument('-t', '--to', type=str, required=False, help='Storage directory')
parser.add_argument('-z', '--zip', type=str, required=False, help='Name of the zip file (including directory if '
                                                                  'necessary)')
parser.add_argument('--delete', action='store_true', required=False, default=False,
                    help='Do you want to delete the folders once transfer is complete?')
args = parser.parse_args()
start_dir = Path(args.start)
args.to=args.start if not args.to else args.to
destination = Path(args.to)
destination.mkdir(parents=True,exist_ok=True)
for item in start_dir.iterdir():
    if item.is_file():
        copy(item, destination)
    if item.is_dir():
        match item.stem:
            case x if x.startswith(ChimeraClasses.ALPHAFOLD_FOLDERNAME):
                [skeletonize_alphafold_folder(folder.__str__(), destination.joinpath(item.stem, folder.stem).__str__()) for folder in
                 item.iterdir()]
            case x if x.startswith(ChimeraClasses.GMX_FOLDERNAME):
                [skeletionize_gmx_folder(folder, destination.joinpath(item.stem, folder.stem)) for folder in
                    item.iterdir() if folder.is_dir()]
            case x if x.startswith(ChimeraClasses.FASTA_FOLDERNAME):
                copytree(item,destination.joinpath(item.stem),dirs_exist_ok=True)
print("Transfer Complete")
if args.zip:
    system(f'zip -mrj {Path(args.zip).with_suffix(".zip")} ')
