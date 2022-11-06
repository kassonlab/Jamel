ProteinList=open('List','r').readlines()
from numpy import empty
from RBDFinder import AlignmentFinder
import AccessiontoFasta
import ChimeraGenerator
import os
from concurrent.futures import ProcessPoolExecutor
#S1= 0-539 AYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLSEFRVYSSANNCTFEYVSQPFLKNLREFVFKNIDGYFKIYSKHTPPQGFSALEPLVDLPIGINITRFQTLLAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSYNYLYRNLKPFERDISTEIYNCYFPLQSYGFQPTVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICA
#RBD= 223-424 TSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSYNYLYRNLKPFERDISTEIYNCYFPLQSYGFQPTVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVN

SequenceofInterest='TSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSYNYLYRNLKPFERDISTEIYNCYFPLQSYGFQPTVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVN'
CommandStart='/gpfs/gpfs0/scratch/jws6pq/BridBridCMfiles/MultimerAlphaFold.sh '
CommandEnd=' /scratch/jws6pq/Notebook/Finished\n'
ProteinsPerSlurm=2
ChimeraOnly='Yes'
Domain="RBD"
FastasforRun=empty(len(ProteinList)*2, dtype=object)
AccessionNumber=[x.split()[0] for x in ProteinList]
ProteinList=[x.split()[-1] for x in ProteinList]
SequenceofInterest=[SequenceofInterest for x in ProteinList]
FirstPartner=['/scratch/jws6pq/BridCMfiles/SARS2.fasta' for x in ProteinList]
SecondPartner=["/gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/"+x+".fasta" for x in ProteinList]
Fastafilenames=[["/gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/"+x+".fasta","/gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/3merSARS2w"+x+Domain+".fasta"] for x in ProteinList]
Alignmentfile=['/scratch/jws6pq/Notebook/Alignment/'+x+'onSARS2.aln' for x in ProteinList]
SpliceBoundary1,SpliceBoundary2=[223 for x in ProteinList],[424 for x in ProteinList]
with ProcessPoolExecutor() as exe:
    exe.map(AccessiontoFasta.FastatoAlignmentFinder(ProteinList))
    SpliceBoundary=exe.map(AlignmentFinder,SequenceofInterest)
    SpliceBoundary3,SpliceBoundary4=[x[0] for x in SpliceBoundary],[x[1] for x in SpliceBoundary]
    FirstPartnerSequence=exe.map(ChimeraGenerator.SequenceSplice,FirstPartner,SpliceBoundary1,SpliceBoundary2)
    SecondPartnerSequence=exe.map(ChimeraGenerator.SequenceSplice,FirstPartner,SpliceBoundary3,SpliceBoundary4)
    if len(FirstPartnerSequence[0])==3:
        FirstPartnerSequences=[x[-1] for x in FirstPartnerSequence]
    elif len(FirstPartnerSequence[0])==4:
        FirstPartnerSequences=[[x[2],x[3]] for x in FirstPartnerSequence]
    if len(SecondPartnerSequence[0])==3:
        SecondPartnerSequences=[x[-1] for x in SecondPartnerSequence]
    elif len(SecondPartnerSequence[0])==4:
        SecondPartnerSequences=[[x[2],x[3]] for x in SecondPartnerSequence]

    #iterate through list in steps?
k=0
for line in ProteinList:
    ProteinInfo=line.split()
    AccessiontoFasta.AccessionNumbertoFasta(ProteinInfo[-1],ProteinInfo[0])
    AccessiontoFasta.FastatoAlignmentFinder(ProteinInfo[-1])
    SpliceBoundary=AlignmentFinder('/scratch/jws6pq/Notebook/Alignment/'+ProteinInfo[-1]+'onSARS2.aln',SequenceofInterest='AYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLSEFRVYSSANNCTFEYVSQPFLKNLREFVFKNIDGYFKIYSKHTPPQGFSALEPLVDLPIGINITRFQTLLAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSYNYLYRNLKPFERDISTEIYNCYFPLQSYGFQPTVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICA')
    FastasforRun[k]=ChimeraGenerator.('SARS2.fasta',ProteinInfo[-1]+'.fasta',SpliceBoundary[0],SpliceBoundary[1],Boundary1=1,Boundary2=540,Domain='S1')[0]
    FastasforRun[k+1]=ChimeraGenerator.DomainExchange('SARS2.fasta',ProteinInfo[-1]+'.fasta',SpliceBoundary[0],SpliceBoundary[1],Boundary1=1,Boundary2=540,Domain='S1')[1]
    k+=2

Proteinindex=0
Slurmfilenumber=1
while Proteinindex in range(len(FastasforRun)-1):
    os.system('cp /scratch/jws6pq/BridCMfiles/MultimerAlphaFold.slurm /scratch/jws6pq/BridCMfiles/'+str(Slurmfilenumber)+'MultimerAlphaFold.slurm')
    Slurmfile=open('/scratch/jws6pq/BridCMfiles/'+str(Slurmfilenumber)+'MultimerAlphaFold.slurm','a')
    FullCommand=CommandStart
    Slurmfile.write('\n#SBATCH -o /scratch/jws6pq/BridCMfiles/'+str(Slurmfilenumber)+'multimerslurm.out\n#SBATCH -e /scratch/jws6pq/BridCMfiles/'+str(Slurmfilenumber)+'multimerslurm.err\n#Run program\n')
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
    # os.system('sbatch /scratch/jws6pq/BridCMfiles/'+str(Slurmfilenumber)+'MultimerAlphaFold.slurm')
    Slurmfilenumber+=1