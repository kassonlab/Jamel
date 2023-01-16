def generate_alphafold_files(FolderDir, AlphafoldFolder, PDBDest_Name='NA', PlddtDest_Name='NA'):
    from pickle import load as pload
    from  json import load as jload
    from os import path
    from shutil import copy
    from numpy import savetxt
    Fullpath=FolderDir+AlphafoldFolder+'/'
    if path.exists(Fullpath+'ranking_debug.json'):
        if PDBDest_Name!='NA':
            copy(Fullpath+'ranked_0.pdb',PDBDest_Name)
        if PlddtDest_Name!='NA':
            with open(Fullpath+'ranking_debug.json', 'r') as jfile:
                HighestRankModel=jload(jfile)['order'][0][0:7]
                with open(f'{Fullpath}result_{HighestRankModel}_multimer_v2_pred_0.pkl', 'rb') as pkl:
                    data=pload(pkl)
                    savetxt(PlddtDest_Name,data['plddt'],fmt="%s",delimiter=" ")


def limited_alphafold_transfer(folder_direc,alphafold_folder,pdb_destin,plddt_destin):
    from pickle import load as pload
    from os import path,listdir
    from shutil import copy
    from numpy import save
    full_path=folder_direc+alphafold_folder+'/'
    if path.exists(full_path+'ranking_debug.json'):
        copy(full_path+'ranking_debug.json',plddt_destin+alphafold_folder+'_ranking_debug.json')
        for file in [file for file in listdir(full_path) if file.startswith('ranked')]:
            copy(full_path+file,pdb_destin+alphafold_folder+'_'+file)
        for index,file in enumerate([file for file in listdir(full_path) if file.startswith('result')]):
            with open(full_path+file, 'rb') as pkl:
                data=pload(pkl)
                print(f'{plddt_destin}{alphafold_folder}_{index}_plddt.npy')
                save(f'{plddt_destin}{alphafold_folder}_{index}_plddt.npy', data['plddt'])