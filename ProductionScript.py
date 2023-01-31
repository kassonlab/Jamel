import ChimeraGenerator
import AccessiontoAlignment
from sys import argv
from json import load
""" This script is made for the automated creation of chimeric proteins between one designated protein called a reference protein,
and a designated section within the reference protein sequence, called sequence_of_interest being replaced by homologous sequences outlined in a list
accession number sequences in a protein_info file"""
# 3 command line inputs are required for this production script, a .tsv file with accession numbers in column one
# and the corresponding name you want associated with that sequence, these name will be used throughout all naming
# including:fasta file names, alphafold out folders, in file fasta identifiers, plddt data files and so on
# DO NOT INCLUDE YOUR REFERENCE PROTEIN IN THIS FILE
protein_info=str(argv[1])
# The second is a chimera_arguments.json file modified to your liking, this file has naming conventions for all important
# files generated like fasta and msas, for analysis later keep naming conventions consistent, the naming works by replacing the
# character_to_replace which is defaulted with an asterisk ('*'), with the names supplied in the protein_info list that was generated
# for the last input, there are also a couple operation that come as arrays with enclosed in [] and the first input is just
# quotations "", if you place a # within the quotes it will prevent that file from being created, for example you prevent
# the msa from being calculated by the script using this method like this ["#","/gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/CoronavirusMSA.aln"]
# PROVIDE YOUR REFERENCE IN THE JSON INPUT FILE
argument_json=str(argv[2])
# Lastly and most simply, create a fasta file with the sequence that you want spliced out and replaced of your reference protein
sequence_of_interest_fasta=argv[3]
# These lines are extracting all the inputs and info outlined in the previous comments
with open(sequence_of_interest_fasta, 'r') as fasta:
    sequence_of_interest=''.join([x for x in fasta if x[0] != '>' if x != '']).strip().replace('\n', '')
with open(argument_json, 'rb') as jfile:
    argument_dict=load(jfile)["arguments"]
# In your protein_info file,
with open(protein_info, 'r') as info_list:
    info_list=info_list.readlines()
# column 1: accession numbers
    accession_number = [x.split()[0] for x in info_list]
# column 2: naming conventions per accession number sequence
    protein_list = [x.split()[-1] for x in info_list]
# Here, important information from your command line inputs are turned into list iterables to be used in the map function later
sequence_of_interest=[sequence_of_interest for x in protein_list]
reference_protein=[argument_dict['reference_protein'] for x in protein_list]
msa_file=[argument_dict['msa_file_name'][1] for x in protein_list]
character_to_replace=argument_dict['character_to_replace']
subunits=[argument_dict['number_of_subunits'] for x in protein_list]
email=[argument_dict['email_for_accession'] for x in protein_list]
monomer_fastas=[argument_dict['monomer_fasta'].replace(character_to_replace, protein) for protein in protein_list]
multimer_fastas=[argument_dict['multimer_fasta'].replace(character_to_replace, protein) for protein in protein_list]
chimera_fastas=[argument_dict['chimera_fastas'].replace(character_to_replace,protein) for protein in protein_list]
msa_fasta=argument_dict['msa_fasta']
reference_protein_name=[argument_dict['reference_protein_fasta_identifier'] for protein in protein_list]
# Here it is being determined whether you require multimer files for your protein, monomers files are created by default,
# matter what for alignment purposes
if subunits==1:
    list(map(AccessiontoAlignment.accession_to_fasta, monomer_fastas, accession_number, email, subunits))
else:
    list(map(AccessiontoAlignment.accession_to_fasta, monomer_fastas, accession_number, email, subunits, multimer_fastas))
# This is the code enacting the msa creation switch that's called by putting a # in the json file, outlined in above comments
if argument_dict['msa_file_name'][0]!='#':
    AccessiontoAlignment.multiple_sequence_alignment(monomer_fastas, msa_fasta,
                                                     msa_file[0],
                                                     reference_protein[0])
# Here splicing information like the boundaries containing the sequence outlined in sequence_of_interest_fasta for your reference
# protein and the sequence from your partner proteins that aligns with the sequence_of_interest in the msa
splice_info=list(map(AccessiontoAlignment.alignment_finder,msa_file, sequence_of_interest,protein_list,reference_protein_name))
spliced_comparison_sequence=[x[0] for x in splice_info]
reference_splice_boundaries =[x[2] for x in splice_info]
# The part of the reference sequence that remains after splicing out the sequence_of_interest is created and marked for splicing
# with the homologous sequence section
reference_sequence_included=[x[1] for x in map(ChimeraGenerator.sequence_splice, reference_protein, reference_splice_boundaries)]
# The 2 sections from the reference and splice partner are combined and put into a fasta file
chimera_sequences=list(map(ChimeraGenerator.chimera_sequence_creation, spliced_comparison_sequence, reference_sequence_included))
list(map(ChimeraGenerator.fasta_creation, chimera_fastas, chimera_sequences, subunits))
# This is an optional file that can be created that holds a list of all the fasta files that were created, chimeric or otherwise,
# this list can be useful when inputting fasta files to alphafold, especially if you're scripting slurm jobs, this file's
# creation can be turned on by removing the # from within the quotes in the json file on the dictionary key "fasta_file_list_name"
if argument_dict['fasta_file_list_name'][0]!='#':
    if subunits==1: fasta_list= monomer_fastas + chimera_fastas
    else: fasta_list= multimer_fastas + chimera_fastas
    with open(argument_dict['fasta_file_list_name'][1], 'w') as fasta_list_file:
        for fasta in fasta_list:
            fasta_list_file.write(f'{fasta}\n')
