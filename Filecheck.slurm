#!/bin/python3
 
#SBATCH -p standard          # partition
#SBATCH -N 1           # number of nodes
#SBATCH -c 4           # number of cores
#SBATCH --mail-user=jws6pq@virginia.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mem=200000
#SBATCH -t 10:00:00     # time
#SBATCH -o /scratch/jws6pq/CMfiles/filecheck.out
#SBATCH -e /scratch/jws6pq/CMfiles/filecheck.err
#SBATCH -A kas_dev
#Run program

import os
NativeList={'/scratch/jws6pq/Notebook/PDB/3mer'+x.split()[-1]+'.pdb':'/scratch/jws6pq/Notebook/Finished/Fastas/3mer'+x.split()[-1]+'.fasta' for x in open('/scratch/jws6pq/CMfiles/List', 'r').readlines()}
ChimeraList={'/scratch/jws6pq/Notebook/PDB/3merSARS2w'+x.split()[-1]+'S1.pdb':'/scratch/jws6pq/Notebook/Finished/Fastas/3merSARS2w'+x.split()[-1]+'S1.fasta' for x in open('/scratch/jws6pq/CMfiles/List', 'r').readlines()}
Fulllist={**ChimeraList,**NativeList}
Submitlist=[]
# for file in Fulllist.keys():
#     if os.path.exists(file)==False:
#         Submitlist.append(Fulllist[file])
Submitlist=[Fulllist[file] for file in Fulllist.keys() if os.path.exists(file)==False]
CommandStart='/scratch/jws6pq/CMfiles/MultimerAlphaFold.sh '
CommandEnd=' /scratch/jws6pq/Notebook/Finished'
# ProteinsPerSlurm=3
# k=0
# Proteinindex=0
SubmitList=['/scratch/jws6pq/BridNotebook/Fastas/3mer'+x.split()[-1]+'.fasta' for x in open('/scratch/jws6pq/BridCMfiles/List').readlines()]
Submit=[SubmitList.append('/scratch/jws6pq/BridNotebook/Fastas/3merSARS2w'+x.split()[-1]+'S1.fasta') for x in open('/scratch/jws6pq/BridCMfiles/List').readlines()]
Slurmfilenumber=0

for file in SubmitList:
    os.system('cp /scratch/jws6pq/BridCMfiles/MultimerAlphaFold.slurm /scratch/jws6pq/BridCMfiles/'+str(Slurmfilenumber)+'MultimerAlphaFold.slurm')
    Slurmfile=open('/scratch/jws6pq/BridCMfiles/'+str(Slurmfilenumber)+'MultimerAlphaFold.slurm','a')
    FullCommand=CommandStart+file+CommandEnd
    Slurmfile.write('\n#SBATCH -o /scratch/jws6pq/BridCMfiles/'+str(Slurmfilenumber)+'multimerslurm.out\n#Run program\n'+FullCommand)
    Slurmfile.close()
    os.system('sbatch /scratch/jws6pq/BridCMfiles/'+str(Slurmfilenumber)+'MultimerAlphaFold.slurm')
    Slurmfilenumber+=1
