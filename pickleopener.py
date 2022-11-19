def GeneratePDBorPlddtfromAlphafold(FolderDir,AlphafoldFolder,PDBDest_Name='NA',PlddtDest_Name='NA'):
    from pickle import load as pload
    from  json import load as jload
    from os import path
    from shutil import copy
    from numpy import savetxt
    Fullpath=FolderDir+AlphafoldFolder+'/'
    if path.exists(Fullpath+'ranking_debug.json'):
        if PDBDest_Name!='NA' and path.exists(PDBDest_Name)==False:
            copy(Fullpath+'ranked_0.pdb',PDBDest_Name)
        if PlddtDest_Name!='NA' and path.exists(PlddtDest_Name)==False:
            jfile = open(Fullpath+'ranking_debug.json', 'r')
            data = jload(jfile)
            HighestRankModel = data['order'][0][0:7]
            infile = open(Fullpath+'result_' + HighestRankModel + '_multimer_v2_pred_0.pkl', 'rb')
            dat = pload(infile)
            savetxt(PlddtDest_Name,dat['plddt'],fmt="%s",delimiter=" ")