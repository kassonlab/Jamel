from Chimeragenesis import ChimeraGenerator
from os import system
SpliceLength1=200
# Protein1=['Chimera_sequeincing.fasta' for x in range(15)]
Boundary1=[x for x in range(0,489-200,50)]
Boundary2=[x+SpliceLength1 for x in Boundary1]
SequenceList1=list(map(ChimeraGenerator.sequence_splice, Protein1, Boundary1, Boundary2))
SpliceLength1=200
Protein2=['SARS2.fasta' for x in range(15)]
Boundary3=[x for x in range(540,973-200,50)]
Boundary4=[x+SpliceLength1 for x in Boundary3]
SequenceList2=list(map(ChimeraGenerator.sequence_splice, Protein2, Boundary3, Boundary4))
Slurmfilenumber=1
for i in range(len(Boundary1)):
    for j in range(len(Boundary3)):
        ChimeraSequence=SequenceList1[i][0]+SequenceList2[j][0]
        Filename='/gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/HA' + str(Boundary1[i]+1) + 'to' + str(Boundary2[i]+1) + 'Spike' + str(Boundary3[j]+1) + 'to' + str(Boundary4[j]+1) + '.fasta'
        ChimeraGenerator.fasta_creation(Filename, ChimeraSequence, 3)
        system('cp /scratch/jws6pq/BridCMfiles/MultimerAlphaFold.slurm /scratch/jws6pq/BridCMfiles/' + str(Slurmfilenumber) + 'HASPIKEAlphaFold.slurm')
        Slurmfile = open('/scratch/jws6pq/BridCMfiles/' + str(Slurmfilenumber) + 'HASPIKEAlphaFold.slurm', 'a')
        Slurmfile.write('\n#SBATCH -e /scratch/jws6pq/BridCMfiles/' + str(Slurmfilenumber) + 'multimerslurm.out\n#Run program\n')
        Slurmfile.write('/gpfs/gpfs0/scratch/jws6pq/BridCMfiles/MultimerAlphaFold.sh '+Filename+' /scratch/jws6pq/Notebook/Finished')
        Slurmfile.close()
        # system('sbatch /scratch/jws6pq/BridCMfiles/' + str(Slurmfilenumber) + 'HASPIKEAlphaFold.slurm')
        Slurmfilenumber+=1

