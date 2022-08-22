
import pickle
import numpy as np
import json
jfile=open('ranking_debug.json','r')
data=json.load(jfile)
HighestRankModel=data['order'][0][0:7]
infile=open('result_'+HighestRankModel+'_multimer_v2_pred_0.pkl','rb')
dat=pickle.load(infile)
np.savetxt('plddt.txt',dat['plddt'],fmt="%s",delimiter=" ")
