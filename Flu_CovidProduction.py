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
        print(f'HA{Boundary1[i] + 1}to{Boundary2[i]+1}Spike{Boundary3[j] + 1}to{Boundary4[j]+1}')
        # Filename=f'/gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/HA{Boundary1[i] + 1}to{Boundary2[i]+1}Spike{Boundary3[j] + 1}to{Boundary4[j]+1}.fasta'
        native_fasta_one=f'/gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/HA{Boundary1[i] + 1}to{Boundary2[i]+1}.fasta'
        native_fasta_two=f'/gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/Spike{Boundary3[j] + 1}to{Boundary4[j]+1}.fasta'
        # ChimeraGenerator.fasta_creation(Filename, ChimeraSequence)
        ChimeraGenerator.fasta_creation(native_fasta_one,SequenceList1[i])
        ChimeraGenerator.fasta_creation(native_fasta_two,SequenceList2[j])
        # Filenames.append(Filename,native_fasta_one,native_fasta_two)
        Filenames.extend([native_fasta_one, native_fasta_two])
Filenames=list(set(Filenames))
for splice1,splice2 in zip(Boundary1,Boundary2):
    print(f'HA{splice1 + 1}to{splice2+1}')
for splice1,splice2 in zip(Boundary3,Boundary4):
    print(f'Spike{splice1 + 1}to{splice2+1}')
Proteinsperslurm=15
fileindex=0
while fileindex in range(len(Filenames)):
    system(f'cp /gpfs/gpfs0/scratch/jws6pq/BridCMfiles/Flu_CovidAlphaFold.slurm /scratch/jws6pq/BridCMfiles/{Slurmfilenumber}Flu_CovidAlphaFold.slurm')
    Slurmfile = open(f'/scratch/jws6pq/BridCMfiles/{Slurmfilenumber}Flu_CovidAlphaFold.slurm', 'a')
    Slurmfile.write(f'\n#SBATCH -e /scratch/jws6pq/BridCMfiles/{Slurmfilenumber}Flu_Covidslurm.out\n#Run program\n')
    proteins_to_run=','.join(Filenames[fileindex:fileindex+Proteinsperslurm])
    Slurmfile.write(f'/gpfs/gpfs0/scratch/jws6pq/BridCMfiles/MultimerAlphaFold.sh {proteins_to_run} /scratch/jws6pq/Notebook/AlphaFold')
    Slurmfile.close()
    system(f'sbatch /scratch/jws6pq/BridCMfiles/{Slurmfilenumber}Flu_CovidAlphaFold.slurm')
    fileindex+=Proteinsperslurm
    Slurmfilenumber+=1

