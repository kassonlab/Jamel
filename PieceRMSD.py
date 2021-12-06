import numpy as np

file = open("/sfs/lustre/bahamut/scratch/jws6pq/Notebook/Finished/Fastas/List", "r")
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
cmd.load("/sfs/lustre/bahamut/scratch/jws6pq/CMfiles/SARS2.pdb")

cmd.remove('organic')
for line in List:
    x = line.split()
    protein = '/sfs/lustre/bahamut/scratch/jws6pq/Notebook/Finished/'+x[-1]+'/' + x[-1] + '.pdb'
    cmd.load(protein)
    i += 1

#You have to manually set the domains, every spike is slightly different
cmd.select('SARS2NTD', selection='SARS2 and resi 1-233')
cmd.select('SARS2RBD', selection='SARS2 and resi 234-432')
cmd.select('SARS2Stalk', selection='SARS2 and resi 433-973')
i=0
for line in List:
    x = line.split()
    RBD=NTD=Stalk =x[-1]
    RMSDs[i+1,0]=x[-1]
    NTD += ' and resi ' + str(number) +'-233'
    RBD += ' and resi 234-432'
    Stalk += ' and resi 433-973'
    NTDName = x[-1] +'ntd'
    RBDName = x[-1] +'rbd'
    StalkName = x[-1] +'stalk'
    cmd.select(NTDName, selection=NTD)
    RMSD=cmd.align('SARS2NTD', x[-1])
    RMSDs[i+1, 1] = RMSD[0]
    cmd.select(RBDName, selection=RBD)
    RMSD = cmd.align('SARS2RBD', x[-1])
    RMSDs[i+1, 2] = RMSD[0]
    cmd.select(StalkName, selection=Stalk)
    RMSD = cmd.align('SARS2Stalk', x[-1])
    RMSDs[i+1, 3] = RMSD[0]
    RMSD = cmd.align('FullSARS2', x[-1])
    RMSDs[i+1, 4] = RMSD[0]
    i += 1

np.savetxt('/sfs/lustre/bahamut/scratch/jws6pq/CMfiles/PieceRMSD.tsv', RMSDs, fmt="%s", delimiter=" ")
