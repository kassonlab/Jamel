ProteinList=open('FullList','r').readlines()
from numpy import empty
from RBDFinder import AlignmentFinder
import AccessiontoFasta
import ChimeraGenerator
import os
from concurrent.futures import ProcessPoolExecutor

#S1= these boundaries along with rbdfinder have to be fixed to includde the last residue 0-539 AYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLSEFRVYSSANNCTFEYVSQPFLKNLREFVFKNIDGYFKIYSKHTPPQGFSALEPLVDLPIGINITRFQTLLAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSYNYLYRNLKPFERDISTEIYNCYFPLQSYGFQPTVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICA
#RBD= 223-424 TSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSYNYLYRNLKPFERDISTEIYNCYFPLQSYGFQPTVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVN
#This production script is for making SARS2 chimera. 
#The Second partner's outlined sequence is spliced into the first

SequenceofInterest='TSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSYNYLYRNLKPFERDISTEIYNCYFPLQSYGFQPTVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVN'
Executable,DestinationFolder='/gpfs/gpfs0/scratch/jws6pq/BridCMfiles/MultimerAlphaFold.sh ',' /scratch/jws6pq/Notebook/Finished\n'
ProteinsPerSlurm=2
ChimeraOnly='Yes'
DomainSetting=['RBD' for x in ProteinList]
AccessionNumber=[x.split()[0] for x in ProteinList]
ProteinList=[x.split()[-1] for x in ProteinList]
SequenceofInterest=[SequenceofInterest for x in ProteinList]
FirstPartner=['/scratch/jws6pq/BridCMfiles/SARS2.fasta' for x in ProteinList]
SecondPartner=["/gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/"+x+".fasta" for x in ProteinList]
Fastafilenames=["/gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/3mer"+x+".fasta" for x in ProteinList]
Chimerafilenames=["/gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/3merSARS2w"+x+DomainSetting[0]+".fasta" for x in ProteinList]
Alignmentfile=['/scratch/jws6pq/Notebook/Alignment/'+x+'onSARS2.aln' for x in ProteinList]
Subunits=[3 for x in ProteinList]
#FP is FirstPartner and SP is SecondPartner
FPBoundary1,FPBoundary2=[223 for x in ProteinList],[424 for x in ProteinList]

with ProcessPoolExecutor() as exe:
    # exe.map(AccessiontoFasta.FastatoAlignmentFinder,ProteinList)
    SpliceBoundary=list(exe.map(AlignmentFinder,Alignmentfile,SequenceofInterest))
    SPBoundary3,SPBoundary4=[x[0] for x in SpliceBoundary],[x[1] for x in SpliceBoundary]
    FirstPartnerSequence=list(exe.map(ChimeraGenerator.SequenceSplice,FirstPartner,FPBoundary1,FPBoundary2))
    FirstPartnerSequence=[x[1] for x in FirstPartnerSequence]
    SecondPartnerSequence=list(exe.map(ChimeraGenerator.SequenceSplice,SecondPartner,SPBoundary3,SPBoundary4))
    SecondPartnerSequence=[x[0] for x in SecondPartnerSequence]
    ChimeraSequences=list(exe.map(ChimeraGenerator.ChimeraSequenceCreation,SecondPartnerSequence,FirstPartnerSequence))
    exe.map(ChimeraGenerator.FastaCreation,Chimerafilenames,ChimeraSequences,Subunits)
Fileindex=0
Slurmfilenumber=1
if ChimeraOnly=='Yes':
    while Fileindex in range(len(Chimerafilenames)):
        os.system('cp /scratch/jws6pq/BridCMfiles/MultimerAlphaFold.slurm /scratch/jws6pq/BridCMfiles/' + str(Slurmfilenumber) + 'MultimerAlphaFold.slurm')
        Slurmfile = open('/scratch/jws6pq/BridCMfiles/' + str(Slurmfilenumber) + 'MultimerAlphaFold.slurm', 'a')
        Slurmfile.write('\n#SBATCH -e /scratch/jws6pq/BridCMfiles/' + str(Slurmfilenumber) + 'multimerslurm.out\n#Run program\n')
        files=','.join(Chimerafilenames[Fileindex:Fileindex+ProteinsPerSlurm])
        Slurmfile.write(Executable+files+DestinationFolder)
        Slurmfile.close()
        os.system('sbatch /scratch/jws6pq/BridCMfiles/'+str(Slurmfilenumber)+'MultimerAlphaFold.slurm')
        Fileindex+=ProteinsPerSlurm
        Slurmfilenumber+=1

if ChimeraOnly=='No':
    Fullfilelist=Fastafilenames+Chimerafilenames
    while Fileindex in range(len(Fullfilelist)):
        os.system('cp /scratch/jws6pq/BridCMfiles/MultimerAlphaFold.slurm /scratch/jws6pq/BridCMfiles/' + str(Slurmfilenumber) + 'MultimerAlphaFold.slurm')
        Slurmfile = open('/scratch/jws6pq/BridCMfiles/' + str(Slurmfilenumber) + 'MultimerAlphaFold.slurm', 'a')
        Slurmfile.write('\n#SBATCH -e /scratch/jws6pq/BridCMfiles/' + str(Slurmfilenumber) + 'multimerslurm.out\n#Run program\n')
        files=','.join(Chimerafilenames[Fileindex:Fileindex+ProteinsPerSlurm])
        Slurmfile.write(Executable+files+DestinationFolder)
        Slurmfile.close()
        os.system('sbatch /scratch/jws6pq/BridCMfiles/'+str(Slurmfilenumber)+'MultimerAlphaFold.slurm')
        Fileindex+=ProteinsPerSlurm
        Slurmfilenumber+=1