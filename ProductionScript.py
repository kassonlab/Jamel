
import ChimeraGenerator
# from os import system as syst
from concurrent.futures import ProcessPoolExecutor
import AccessiontoAlignment
#S1= these boundaries along with rbdfinder have to be fixed to includde the last residue 0-539 AYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLSEFRVYSSANNCTFEYVSQPFLKNLREFVFKNIDGYFKIYSKHTPPQGFSALEPLVDLPIGINITRFQTLLAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSYNYLYRNLKPFERDISTEIYNCYFPLQSYGFQPTVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICA
#This production script is for making S1 SARS2 chimera.
#The Second partner's outlined sequence is spliced into the first
info_list=open("/gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/List_of_Coronaviruses",'r').readlines()
sequence_of_interest= 'AYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDP' \
                   'FLSEFRVYSSANNCTFEYVSQPFLKNLREFVFKNIDGYFKIYSKHTPPQGFSALEPLVDLPIGINITRFQTLLAAYYVGYLQPRTFLLKYNENGTI' \
                   'TDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYG' \
                   'VSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSYNYLYRNLKPFERDISTEIYNCYFPLQSYGFQPT' \
                   'VGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVI' \
                   'TPGTNTSNQVAVLYQDVNCTEVTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICA'
executable, destination_folder= '/gpfs/gpfs0/scratch/jws6pq/BridCMfiles/MultimerAlphaFold.sh ',\
                             ' /scratch/jws6pq/Notebook/Finished\n'
proteins_per_slurm=2
chimera_only= 'Yes'
domain_setting=['S1' for protein in info_list]
accession_number=[x.split()[0] for x in info_list]
protein_list=[x.split()[-1] for x in info_list]
sequence_of_interest=[sequence_of_interest for x in protein_list]
SARS2=['/scratch/jws6pq/BridCMfiles/SARS2.fasta' for x in protein_list]
monomer_fasta=[f"/gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/{protein}.fasta" for protein in protein_list]
multimer_fasta=[f"/gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/3mer{protein}.fasta" for protein in protein_list]
chimera_fastas=[f"/gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/3merSARS2w{protein}{domain_setting[0]}.fasta" for protein in protein_list]
subunits=[3 for x in protein_list]
email=['jws6pq@virginia.edu' for x in protein_list]
msa_file=['/gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/CoronavirusMSA.aln' for x in protein_list]
# SP is splice_partner
sars_boundary_one, sars_boundary_two= [0 for x in protein_list], [540 for x in protein_list]

with ProcessPoolExecutor(max_workers=6) as exe:
    exe.map(AccessiontoAlignment.accession_to_fasta, protein_list, monomer_fasta, accession_number,email,subunits,multimer_fasta)
    AccessiontoAlignment.multiple_sequence_alignment(monomer_fasta,'/gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/CoronavirusMSA.fasta',
                                                     '/gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/CoronavirusMSA.aln')
    sp_sequence=[x[0] for x in exe.map(AccessiontoAlignment.alignment_finder,msa_file, sequence_of_interest,protein_list)]
    sars_sequence=[x[1] for x in exe.map(ChimeraGenerator.sequence_splice, SARS2, sars_boundary_one, sars_boundary_two)]
    chimera_sequences=list(exe.map(ChimeraGenerator.chimera_sequence_creation, sp_sequence, sars_sequence))
    exe.map(ChimeraGenerator.fasta_creation, chimera_fastas, chimera_sequences, subunits)
# Fileindex=0
# Slurmfilenumber=1
# if chimera_only== 'Yes':
#     while Fileindex in range(len(chimera_fastas)):
#         syst(f'cp /scratch/jws6pq/BridCMfiles/MultimerAlphaFold.slurm /scratch/jws6pq/BridCMfiles/{Slurmfilenumber}MultimerAlphaFold.slurm')
#         Slurmfile = open(f'/scratch/jws6pq/BridCMfiles/{Slurmfilenumber}MultimerAlphaFold.slurm', 'a')
#         Slurmfile.write(f'\n#SBATCH -e /scratch/jws6pq/BridCMfiles/{Slurmfilenumber}multimerslurm.out\n#Run program\n')
#         files=','.join(chimera_fastas[Fileindex:Fileindex + proteins_per_slurm])
#         Slurmfile.write(executable + files + destination_folder)
#         Slurmfile.close()
#         syst(f'sbatch /scratch/jws6pq/BridCMfiles/{Slurmfilenumber}MultimerAlphaFold.slurm')
#         Fileindex+=proteins_per_slurm
#         Slurmfilenumber+=1
#
# if chimera_only== 'No':
#     Fullfilelist= fasta_files + chimera_fastas
#     while Fileindex in range(len(Fullfilelist)):
#         syst(f'cp /scratch/jws6pq/BridCMfiles/MultimerAlphaFold.slurm /scratch/jws6pq/BridCMfiles/{Slurmfilenumber}MultimerAlphaFold.slurm')
#         Slurmfile = open(f'/scratch/jws6pq/BridCMfiles/{Slurmfilenumber}MultimerAlphaFold.slurm', 'a')
#         Slurmfile.write(f'\n#SBATCH -e /scratch/jws6pq/BridCMfiles/{Slurmfilenumber}multimerslurm.out\n#Run program\n')
#         files=','.join(chimera_fastas[Fileindex:Fileindex + proteins_per_slurm])
#         Slurmfile.write(executable + files + destination_folder)
#         Slurmfile.close()
#         syst(f'sbatch /scratch/jws6pq/BridCMfiles/{Slurmfilenumber}MultimerAlphaFold.slurm')
#         Fileindex+=proteins_per_slurm
#         Slurmfilenumber+=1