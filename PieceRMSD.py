import numpy as np

file = open("E:\ResearchScripts\list", "r")
List = file.readlines()
Listlength = len(List)
RMSDs = np.empty(((Listlength+1),5), dtype=object)
RMSDs[0,0]='Name'
RMSDs[0,1]='NTD'
RMSDs[0,2]='RBD'
RMSDs[0,3]='Stalk'
RMSDs[0,4]='Total'
number=1
i = 0
cmd.load("E:\Research\\FullSARS2.pdb")

cmd.remove('organic')
for line in List:
    x = line.split()
    protein = 'E:\Research\\' + x[0] + '.pdb'
    cmd.load(protein)
    i += 1

#You have to manually set the domains, every spike is slightly different
cmd.select('SARS2NTD', selection='FullSARS2 and resi 1-319')
cmd.select('SARS2RBD', selection='FullSARS2 and resi 319-541')
cmd.select('SARS2Stalk', selection='FullSARS2 and resi 541-1273')
i=0
for line in List:
    x = line.split()
    RBD=NTD=Stalk =x[0]
    RMSDs[i+1,0]=x[0]
    NTD += ' and resi ' + str(number) +'-319'
    RBD += ' and resi 319-541'
    Stalk += ' and resi 541-1273'
    NTDName = x[0] +'ntd'
    RBDName = x[0] +'rbd'
    StalkName = x[0] +'stalk'
    cmd.select(NTDName, selection=NTD)
    RMSD=cmd.align('SARS2NTD', x[0])
    RMSDs[i+1, 1] = RMSD[0]
    cmd.select(RBDName, selection=RBD)
    RMSD = cmd.align('SARS2RBD', x[0])
    RMSDs[i+1, 2] = RMSD[0]
    cmd.select(StalkName, selection=Stalk)
    RMSD = cmd.align('SARS2Stalk', x[0])
    RMSDs[i+1, 3] = RMSD[0]
    RMSD = cmd.align('FullSARS2', x[0])
    RMSDs[i+1, 4] = RMSD[0]
    i += 1

np.savetxt('E:\PythonShit\PANGPieceRMSD.tsv', RMSDs, fmt="%s", delimiter=" ")