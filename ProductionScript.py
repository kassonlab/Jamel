from RBDFinder import alignment_finder
import ChimeraGenerator
from os import system as syst
from concurrent.futures import ProcessPoolExecutor

#S1= these boundaries along with rbdfinder have to be fixed to includde the last residue 0-539 AYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLSEFRVYSSANNCTFEYVSQPFLKNLREFVFKNIDGYFKIYSKHTPPQGFSALEPLVDLPIGINITRFQTLLAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSYNYLYRNLKPFERDISTEIYNCYFPLQSYGFQPTVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICA
#This production script is for making S1 SARS2 chimera.
#The Second partner's outlined sequence is spliced into the first
info_list=open('List','r').readlines()

SequenceofInterest='AYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDP' \
                   'FLSEFRVYSSANNCTFEYVSQPFLKNLREFVFKNIDGYFKIYSKHTPPQGFSALEPLVDLPIGINITRFQTLLAAYYVGYLQPRTFLLKYNENGTI' \
                   'TDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYG' \
                   'VSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSYNYLYRNLKPFERDISTEIYNCYFPLQSYGFQPT' \
                   'VGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVI' \
                   'TPGTNTSNQVAVLYQDVNCTEVTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICA'
Executable,DestinationFolder='/gpfs/gpfs0/scratch/jws6pq/BridCMfiles/MultimerAlphaFold.sh ',\
                             ' /scratch/jws6pq/Notebook/Finished\n'
proteins_per_slurm=2
chimera_only= 'Yes'
DomainSetting=['S1' for protein in info_list]
AccessionNumber=[x.split()[0] for x in info_list]
protein_list=[x.split()[-1] for x in info_list]
SequenceofInterest=[SequenceofInterest for x in protein_list]
SARS2=['/scratch/jws6pq/BridCMfiles/SARS2.fasta' for x in protein_list]
splice_partner=[f"/gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/{protein}.fasta" for protein in protein_list]
Fastafilenames=[f"/gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/3mer{protein}.fasta" for protein in protein_list]
Chimerafilenames=[f"/gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/3merSARS2w{protein}{DomainSetting[0]}.fasta" for protein in protein_list]
Alignmentfile=[f'/scratch/jws6pq/Notebook/Alignment/{protein}onSARS2.aln' for protein in protein_list]
Subunits=[3 for x in protein_list]
# SP is splice_partner
sars_boundary_one, sars_boundary_two= [0 for x in protein_list], [540 for x in protein_list]

with ProcessPoolExecutor() as exe:
    # exe.map(AccessiontoFasta.FastatoAlignmentFinder,ProteinList)
    SpliceBoundary=list(exe.map(alignment_finder, Alignmentfile, SequenceofInterest))
    SPBoundary3,SPBoundary4=[x[0] for x in SpliceBoundary],[x[1] for x in SpliceBoundary]
    sars_sequence=list(exe.map(ChimeraGenerator.sequence_splice, SARS2, sars_boundary_one, sars_boundary_two)[1])
    sp_sequence=list(exe.map(ChimeraGenerator.sequence_splice, splice_partner, SPBoundary3, SPBoundary4)[0])
    ChimeraSequences=list(exe.map(ChimeraGenerator.chimera_sequence_creation, sp_sequence, sars_sequence))
    exe.map(ChimeraGenerator.fasta_creation, Chimerafilenames, ChimeraSequences, Subunits)
Fileindex=0
Slurmfilenumber=1
if chimera_only== 'Yes':
    while Fileindex in range(len(Chimerafilenames)):
        syst(f'cp /scratch/jws6pq/BridCMfiles/MultimerAlphaFold.slurm /scratch/jws6pq/BridCMfiles/{Slurmfilenumber}MultimerAlphaFold.slurm')
        Slurmfile = open(f'/scratch/jws6pq/BridCMfiles/{Slurmfilenumber}MultimerAlphaFold.slurm', 'a')
        Slurmfile.write(f'\n#SBATCH -e /scratch/jws6pq/BridCMfiles/{Slurmfilenumber}multimerslurm.out\n#Run program\n')
        files=','.join(Chimerafilenames[Fileindex:Fileindex + proteins_per_slurm])
        Slurmfile.write(Executable+files+DestinationFolder)
        Slurmfile.close()
        syst(f'sbatch /scratch/jws6pq/BridCMfiles/{Slurmfilenumber}MultimerAlphaFold.slurm')
        Fileindex+=proteins_per_slurm
        Slurmfilenumber+=1

if chimera_only== 'No':
    Fullfilelist=Fastafilenames+Chimerafilenames
    while Fileindex in range(len(Fullfilelist)):
        syst(f'cp /scratch/jws6pq/BridCMfiles/MultimerAlphaFold.slurm /scratch/jws6pq/BridCMfiles/{Slurmfilenumber}MultimerAlphaFold.slurm')
        Slurmfile = open(f'/scratch/jws6pq/BridCMfiles/{Slurmfilenumber}MultimerAlphaFold.slurm', 'a')
        Slurmfile.write(f'\n#SBATCH -e /scratch/jws6pq/BridCMfiles/{Slurmfilenumber}multimerslurm.out\n#Run program\n')
        files=','.join(Chimerafilenames[Fileindex:Fileindex + proteins_per_slurm])
        Slurmfile.write(Executable+files+DestinationFolder)
        Slurmfile.close()
        syst(f'sbatch /scratch/jws6pq/BridCMfiles/{Slurmfilenumber}MultimerAlphaFold.slurm')
        Fileindex+=proteins_per_slurm
        Slurmfilenumber+=1