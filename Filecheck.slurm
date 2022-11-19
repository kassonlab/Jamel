#!/bin/python3
 
#SBATCH -p gpu          # partition
#SBATCH --gres=gpu
#SBATCH -N 1           # number of nodes
#SBATCH --ntasks-per-node=6           # number of cores
#SBATCH --mail-user=jws6pq@virginia.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mem=200000
#SBATCH -t 10:00:00     # time
#SBATCH -o /gpfs/gpfs0/scratch/jws6pq/BridCMfiles/Flusstoberun.out
#SBATCH -e /gpfs/gpfs0/scratch/jws6pq/BridCMfiles/filecheck.err
#SBATCH -A kas_dev
#Run program
from os import path
from os import system
Runstragglers='Yes'
Proteinsperslurm=10
#Multimer or Flu_Covid
Slurmname='Flu_Covid'
Input="/gpfs/gpfs0/scratch/jws6pq/BridCMfiles/Flu_CovidList"
Listtocheck=open(Input,'r').readlines()
Listtocheck=[x.rstrip() for x in Listtocheck]
Listtorun=[]

for file in Listtocheck:
  if path.exists("/gpfs/gpfs0/scratch/jws6pq/Notebook/Plddt/"+file+".plddt")==False:
    Listtorun.append('/gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/'+file+'.fasta')
    print(file.split('SARS2w')[-1].replace('RBD','').replace('S1',''))
def ProductionRun(Listtorun,Slurmname,Proteinsperslurm):
  fileindex=0
  Slurmfilenumber=1
  while fileindex in range(len(Listtorun)):
    system('cp /gpfs/gpfs0/scratch/jws6pq/BridCMfiles/'+Slurmname+'AlphaFold.slurm /scratch/jws6pq/BridCMfiles/' + str(Slurmfilenumber) +Slurmname+'AlphaFold.slurm')
    Slurmfile = open('/scratch/jws6pq/BridCMfiles/' + str(Slurmfilenumber) +Slurmname+'AlphaFold.slurm', 'a')
    Slurmfile.write('\n#SBATCH -e /scratch/jws6pq/BridCMfiles/' + str(Slurmfilenumber) +Slurmname+'slurm.out\n#Run program\n')
    Slurmfile.write('/gpfs/gpfs0/scratch/jws6pq/BridCMfiles/MultimerAlphaFold.sh '+','.join(Listtorun[fileindex:fileindex+Proteinsperslurm])+' /scratch/jws6pq/Notebook/Finished')
    Slurmfile.close()
#    system('sbatch /scratch/jws6pq/BridCMfiles/' + str(Slurmfilenumber) +Slurmname+'AlphaFold.slurm')
    fileindex+=Proteinsperslurm
    Slurmfilenumber+=1
if Runstragglers=='Yes':
  ProductionRun(Listtorun,Slurmname,Proteinsperslurm)
