def GeneratePDBandPlddtfromAlphafold(FolderDir,AlphafoldFolder,PDBDest_Name='NA',PlddtDest_Name='NA'):
    import pickle
    import json
    import os
    import shutil
    import numpy as np
    Fullpath=FolderDir+AlphafoldFolder+'/'
    if os.path.exists(Fullpath+'ranking_debug.json'):
        if PDBDest_Name!='NA' and os.path.exists(PDBDest_Name)==False:
            shutil.copy(Fullpath+'ranked_0.pdb',PDBDest_Name)
        if PlddtDest_Name!='NA' and os.path.exists(PlddtDest_Name)==False:
            jfile = open(Fullpath+'ranking_debug.json', 'r')
            data = json.load(jfile)
            HighestRankModel = data['order'][0][0:7]
            infile = open(Fullpath+'result_' + HighestRankModel + '_multimer_v2_pred_0.pkl', 'rb')
            dat = pickle.load(infile)
            np.savetxt(PlddtDest_Name,dat['plddt'],fmt="%s",delimiter=" ")