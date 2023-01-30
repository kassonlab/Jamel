def generate_alphafold_files(output_folder, new_pdb='NA', new_plddt='NA'):
    from pickle import load as pload
    from  json import load as jload
    from os import path
    from shutil import copy
    from numpy import savetxt
    if path.exists(output_folder + 'ranking_debug.json'):
        if new_pdb!= 'NA':
            copy(output_folder + 'ranked_0.pdb', new_pdb)
        if new_plddt!= 'NA':
            with open(output_folder + 'ranking_debug.json', 'r') as jfile:
                HighestRankModel=jload(jfile)['order'][0]
                with open({HighestRankModel}+'.pkl', 'rb') as pfile:
                    data=pload(pfile)
                    savetxt(new_plddt, data['plddt'], fmt="%s", delimiter=" ")


def limited_alphafold_transfer(folder_direc,alphafold_folder,pdb_destin,plddt_destin):
    from pickle import load as pload
    from os import path,listdir
    from shutil import copy
    from numpy import save
    full_path=folder_direc+alphafold_folder+'/'
    if path.exists(full_path+'ranking_debug.json'):
        copy(full_path+'ranking_debug.json',plddt_destin+alphafold_folder+'ranking_debug.json')
        for file in [file for file in listdir(full_path) if file.startswith('ranked')]:
            copy(full_path+file,pdb_destin+alphafold_folder+'_'+file)
        for index,file in enumerate([file for file in listdir(full_path) if file.startswith('result')]):
            with open(full_path+file, 'rb') as pkl:
                data=pload(pkl)
                save(f'{plddt_destin}{alphafold_folder}_{index}_plddt.npy', data['plddt'])