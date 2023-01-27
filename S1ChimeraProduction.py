import ChimeraGenerator
from concurrent.futures import ProcessPoolExecutor
import AccessiontoAlignment
from sys import argv
# This is an example script for generating fasta files for SARS-CoV-2 Chimeras that replace the S1 of SARS and replace it with
# a number of S1 sequences from other coronaviruses

# Create a tsv with a column 1 containing the nicknames for proteins you want in filenames/fasta designations and column 2 with the accession numbers
# This is tsv will be entered here as the first argument on the command line as protein_info
protein_info=argv[1]
# The next argument is whatever domain you are splicing out of your list of proteins, this is only for naming purposes and
# can be taken out
domain_to_splice=argv[2]
# protein_info is opened and split into its 2 parts of accession numbers and protein names
info_list=open(protein_info,'r').readlines()
accession_number=[x.split()[0] for x in info_list]
protein_list=[x.split()[-1] for x in info_list]
# This is the S1 of SARS-CoV-2 6vsb PDB
# have this be fasta file input
sequence_of_interest= 'AYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDP' \
                   'FLSEFRVYSSANNCTFEYVSQPFLKNLREFVFKNIDGYFKIYSKHTPPQGFSALEPLVDLPIGINITRFQTLLAAYYVGYLQPRTFLLKYNENGTI' \
                   'TDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYG' \
                   'VSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSYNYLYRNLKPFERDISTEIYNCYFPLQSYGFQPT' \
                   'VGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVI' \
                   'TPGTNTSNQVAVLYQDVNCTEVTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICA'
# These are all turned into a list because map() function needs an iterable
domain_setting=[str(domain_to_splice) for protein in info_list]
sequence_of_interest=[sequence_of_interest for x in protein_list]
SARS2=['/scratch/jws6pq/BridCMfiles/SARS2.fasta' for x in protein_list]
msa_file=['/gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/CoronavirusMSA.aln' for x in protein_list]
subunits=[3 for x in protein_list]
email=['example@outlook.com' for x in protein_list]
# Setting naming convention for the SARS2 chimeras
monomer_fasta=[f"/gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/{protein}.fasta" for protein in protein_list]
multimer_fasta=[f"/gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/3mer{protein}.fasta" for protein in protein_list]
chimera_fastas=[f"/gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/3merSARS2w{protein}{domain_setting[0]}.fasta" for protein in protein_list]

with ProcessPoolExecutor(max_workers=6) as exe:
    exe.map(AccessiontoAlignment.accession_to_fasta, monomer_fasta, accession_number,email,subunits,multimer_fasta)
    AccessiontoAlignment.multiple_sequence_alignment(monomer_fasta,'/gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/CoronavirusMSA.fasta',
                                                     '/gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/CoronavirusMSA.aln',
                                                     '/scratch/jws6pq/BridCMfiles/SARS2.fasta')
    splice_info=list(exe.map(AccessiontoAlignment.alignment_finder,msa_file, sequence_of_interest,protein_list,['SARS2' for x in protein_list]))
    sp_sequence=[x[0] for x in splice_info]
    sars_boundary_one =[x[3] for x in splice_info]
    sars_boundary_two=[x[4] for x in splice_info]
    sars_sequence=[x[1] for x in exe.map(ChimeraGenerator.sequence_splice, SARS2, sars_boundary_one, sars_boundary_two)]
    chimera_sequences=list(exe.map(ChimeraGenerator.chimera_sequence_creation, sp_sequence, sars_sequence))
    exe.map(ChimeraGenerator.fasta_creation, chimera_fastas, chimera_sequences, subunits)

