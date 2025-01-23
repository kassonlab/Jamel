from json import load as j_load
from numpy import savetxt
from pickle import load as p_load
from os import path, listdir
from shutil import copy
from numpy import save
def generate_alphafold_files(output_folder, new_plddt='NA', new_pdb='NA'):
    """Creates a text file containing the plddt values of the highest_rank_model extracted from alphafold's result pkl file
    and renames the ranked_0.pdb file and places it in the desired directory."""
    # Checking to see if ranking_debug.json exists. This file is the last to be output by alphafold and is a check that
    # the pkl file you want to extract from exists, as well as to avoid errors
    if path.exists(output_folder + 'ranking_debug.json'):
        if new_pdb != 'NA':
            # The highest ranked structure is copied with a new name and directory
            copy(output_folder + 'ranked_0.pdb', new_pdb)
        if new_plddt != 'NA':
            # ranking_debug is also useful for determining which result pkl file is the highest ranked. The model result pkl files are
            # numbered by the order they are created and not their overall confidence score. The information about their rank by
            # confidence score is found in ranking_debug.json
            with open(output_folder + 'ranking_debug.json', 'r') as jfile:
                highest_rank_model = j_load(jfile)['order'][0]
                with open(f'{output_folder}result_{highest_rank_model}.pkl', 'rb') as pfile:
                    data = p_load(pfile)
                    # The plddt scores are put into a column in a text file named by new_plddt
                    savetxt(new_plddt, data['plddt'], fmt='%s', delimiter=' ')

# def limited_alphafold_transfer(alphafold_dir,storage_dir):
#
#     if path.exists(full_path+'ranking_debug.json'):
#         copy(full_path+'ranking_debug.json',plddt_destin+alphafold_folder+'ranking_debug.json')
#         for file in [file for file in listdir(full_path) if file.startswith('ranked')]:
#             copy(full_path+file,pdb_destin+alphafold_folder+'_'+file)
#         for index,file in enumerate([file for file in listdir(full_path) if file.startswith('result')]):
#             with open(full_path+file, 'rb') as pkl:
#                 data=pload(pkl)
#                 save(f'{plddt_destin}{alphafold_folder}_{index}_plddt.npy', data['plddt'])
target='/gpfs/gpfs0/scratch/jws6pq/Notebook/AlphaFold_Outputs/NTD_Head/3merNTDwA_swine_Oklahoma_A02712158_2022Stalk/'
with open(target + 'ranking_debug.json', 'r') as jfile:
    highest_rank_model = j_load(jfile)['order'][0]
    with open(f'{target}result_{highest_rank_model}.pkl', 'rb') as pfile:
        data = p_load(pfile)
        # The plddt scores are put into a column in a text file named by new_plddt
        savetxt(target+'full.output', data, delimiter=' ')