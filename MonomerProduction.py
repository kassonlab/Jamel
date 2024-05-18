ProteinList=open('List','r').readlines()
import numpy as np
import ChimeraGenerator
import os
CommandStart='python /scratch/jws6pq/CMfiles/multimeralpha_list.py -s='
CommandEnd=' -o=/scratch/jws6pq/Notebook/Finished\n'
ProteinsPerSlurm=6
ChimeraOnly='No'
FastasforRun=np.empty(len(ProteinList)*2, dtype=object)
k=0
for line in ProteinList:
    ProteinInfo=line.split()
    AccessiontoFasta.accession_to_fasta(ProteinInfo[-1], ProteinInfo[0])
    AccessiontoFasta.fasta_to_alignment(ProteinInfo[-1])
    SpliceBoundary=RBDFinder.alignment_finder('/scratch/jws6pq/Notebook/Alignment/' + ProteinInfo[-1] + 'onSARS2.aln', sequence_of_interest='AYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLSEFRVYSSANNCTFEYVSQPFLKNLREFVFKNIDGYFKIYSKHTPPQGFSALEPLVDLPIGINITRFQTLLAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSYNYLYRNLKPFERDISTEIYNCYFPLQSYGFQPTVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICA')
    FastasforRun[k]=ChimeraGenerator.DomainExchange('SARS2.fasta',ProteinInfo[-1]+'.fasta',SpliceBoundary[0],SpliceBoundary[1],Boundary1=1,Boundary2=540,Domain='S1')[0]
    FastasforRun[k+1]=ChimeraGenerator.DomainExchange('SARS2.fasta',ProteinInfo[-1]+'.fasta',SpliceBoundary[0],SpliceBoundary[1],Boundary1=1,Boundary2=540,Domain='S1')[1]
    k+=2
Proteinindex=0
Slurmfilenumber=0
while Proteinindex in range(len(FastasforRun)-1):
    os.system('cp /scratch/jws6pq/CMfiles/MultimerAlphaFold.slurm /scratch/jws6pq/CMfiles/'+str(Slurmfilenumber)+'MultimerAlphaFold.slurm')
    Slurmfile=open('/scratch/jws6pq/CMfiles/'+str(Slurmfilenumber)+'MultimerAlphaFold.slurm','a')
    FullCommand=CommandStart
    Slurmfile.write('#SBATCH -o /scratch/jws6pq/CMfiles/'+str(Slurmfilenumber)+'multimerslurm.out\n#SBATCH -e /scratch/jws6pq/CMfiles/'+str(Slurmfilenumber)+'multimerslurm.err\n#Run program\n')
    Proteinsperslurmindex = 0
    if ChimeraOnly=='No':
        while Proteinsperslurmindex in range(ProteinsPerSlurm):
            FullCommand+=FastasforRun[Proteinindex]
            Proteinindex+=1
            Proteinsperslurmindex+=1
    elif ChimeraOnly=='Yes':
        while Proteinsperslurmindex in range(ProteinsPerSlurm):
            FullCommand += FastasforRun[Proteinindex]
            Proteinindex += 2
            Proteinsperslurmindex += 1
    FullCommand+=CommandEnd
    Slurmfile.write(FullCommand)
    Slurmfile.close()
    # os.system('sbatch /scratch/jws6pq/CMfiles/'+str(Slurmfilenumber)+'MultimerAlphaFold.slurm')
    Slurmfilenumber+=1