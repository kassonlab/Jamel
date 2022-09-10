def GeneratePDBandPlddtfromAlphafold(FolderDir,AlphafoldFolder,PDBDest_Name='NA',PlddtDest_Name='NA'):
    import pickle
    import json
    import os
    import shutil
    import numpy as np
    os.chdir(FolderDir+AlphafoldFolder)
    #os.path.exists(PDBDest_Name)==False and os.path.exists(PlddtDest_Name)==False
    if os.path.exists('ranking_debug.json') and :
        jfile=open('ranking_debug.json','r')
        print(AlphafoldFolder)
        data=json.load(jfile)
        HighestRankModel=data['order'][0][0:7]
        infile=open('result_'+HighestRankModel+'_multimer_v2_pred_0.pkl','rb')
        dat=pickle.load(infile)
        if PlddtDest_Name!='NA':
            np.savetxt(PlddtDest_Name,dat['plddt'],fmt="%s",delimiter=" ")
        if PDBDest_Name!='NA':
            shutil.copy('ranked_0.pdb',PDBDest_Name)

