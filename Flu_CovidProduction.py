#!/bin/python3
#SBATCH -o Flu_CovidList
#SBATCH -p standard
import sys
sys.path.append('/scratch/jws6pq/BridCMfiles')
import ChimeraGenerator
from os import system

SpliceLength1=200
Boundary1=[x for x in range(0,489-SpliceLength1,10)]
Protein1=['HA.fasta' for x in Boundary1]
Boundary2=[x+SpliceLength1 for x in Boundary1]
SequenceList1=[x[0] for x in list(map(ChimeraGenerator.sequence_splice, Protein1, Boundary1, Boundary2))]

SpliceLength2=200
Boundary3=[x for x in range(540,973-SpliceLength2,20)]
Protein2=['SARS2.fasta' for x in Boundary3]
Boundary4=[x+SpliceLength1 for x in Boundary3]
SequenceList2=[x[0] for x in list(map(ChimeraGenerator.sequence_splice, Protein2, Boundary3, Boundary4))]

Filenames=[]
Slurmfilenumber=1
for i in range(len(Boundary1)):
    for j in range(len(Boundary3)):
        ChimeraSequence=SequenceList1[i]+SequenceList2[j]
        print('HA'+str(Boundary1[i]+1) + 'to' + str(Boundary2[i]+1) + 'Spike' + str(Boundary3[j]+1) + 'to' + str(Boundary4[j]+1))
        Filename='/gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/HA' + str(Boundary1[i]+1) + 'to' + str(Boundary2[i]+1) + 'Spike' + str(Boundary3[j]+1) + 'to' + str(Boundary4[j]+1) + '.fasta'
        ChimeraGenerator.fasta_creation(Filename, ChimeraSequence, 3)
        Filenames.append(Filename)
        
Proteinsperslurm=8
fileindex=0
while fileindex in range(len(Filenames)):
    system('cp /gpfs/gpfs0/scratch/jws6pq/BridCMfiles/Flu_CovidAlphaFold.slurm /scratch/jws6pq/BridCMfiles/' + str(Slurmfilenumber) + 'Flu_CovidAlphaFold.slurm')
    Slurmfile = open('/scratch/jws6pq/BridCMfiles/' + str(Slurmfilenumber) + 'Flu_CovidAlphaFold.slurm', 'a')
    Slurmfile.write('\n#SBATCH -e /scratch/jws6pq/BridCMfiles/' + str(Slurmfilenumber) + 'Flu_Covidslurm.out\n#Run program\n')
    Slurmfile.write('/gpfs/gpfs0/scratch/jws6pq/BridCMfiles/MultimerAlphaFold.sh '+','.join(Filenames[fileindex:fileindex+Proteinsperslurm])+' /scratch/jws6pq/Notebook/Finished')
    Slurmfile.close()
    system('sbatch /scratch/jws6pq/BridCMfiles/' + str(Slurmfilenumber) + 'Flu_CovidAlphaFold.slurm')
    fileindex+=Proteinsperslurm
    Slurmfilenumber+=1

