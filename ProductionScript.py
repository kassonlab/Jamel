import ChimeraGenerator
import AccessiontoAlignment
from sys import argv
from json import load
protein_info=str(argv[1])
argument_json=str(argv[2])
sequence_of_interest_fasta=argv[3]

with open(sequence_of_interest_fasta, 'r') as fasta:
    sequence_of_interest=''.join([x for x in fasta if x[0] != '>' if x != '']).strip().replace('\n', '')
# protein_info is opened and split into its 2 parts of accession numbers and protein names
with open(argument_json, 'rb') as jfile:
    argument_dict=load(jfile)["arguments"]
with open(protein_info, 'r') as info_list:
    info_list=info_list.readlines()
    accession_number = [x.split()[0] for x in info_list]
    protein_list = [x.split()[-1] for x in info_list]

sequence_of_interest=[sequence_of_interest for x in protein_list]
reference_protein=[argument_dict['reference_protein'] for x in protein_list]
msa_file=[argument_dict['msa_file_name'][1] for x in protein_list]
character_to_replace=argument_dict['character_to_replace']
subunits=[argument_dict['number_of_subunits'] for x in protein_list]
email=[argument_dict['email_for_accession'] for x in protein_list]
monomer_fasta=[argument_dict['monomer_fasta'].replace(character_to_replace,protein) for protein in protein_list]
multimer_fasta=[argument_dict['multimer_fasta'].replace(character_to_replace,protein) for protein in protein_list]
chimera_fastas=[argument_dict['chimera_fastas'].replace(character_to_replace,protein) for protein in protein_list]
msa_fasta=argument_dict['msa_fasta']
reference_protein_name=[argument_dict['reference_protein_fasta_identifier'] for protein in protein_list]

if subunits==1:
    list(map(AccessiontoAlignment.accession_to_fasta, monomer_fasta, accession_number,email,subunits))
else:
    list(map(AccessiontoAlignment.accession_to_fasta, monomer_fasta, accession_number, email, subunits,multimer_fasta))
if argument_dict['msa_file_name'][0]!='#':
    AccessiontoAlignment.multiple_sequence_alignment(monomer_fasta,msa_fasta,
                                                     msa_file[0],
                                                     reference_protein[0])
splice_info=list(map(AccessiontoAlignment.alignment_finder,msa_file, sequence_of_interest,protein_list,reference_protein_name))
spliced_comparison_sequence=[x[0] for x in splice_info]
reference_splice_boundaries =[x[2] for x in splice_info]
sars_sequence=[x[1] for x in map(ChimeraGenerator.sequence_splice, reference_protein, reference_splice_boundaries)]
chimera_sequences=list(map(ChimeraGenerator.chimera_sequence_creation, spliced_comparison_sequence, sars_sequence))
list(map(ChimeraGenerator.fasta_creation, chimera_fastas, chimera_sequences, subunits))

if argument_dict['fasta_file_list_name'][0]!='#':
    if subunits==1: fasta_list=monomer_fasta+chimera_fastas
    else: fasta_list=multimer_fasta+chimera_fastas
    with open(argument_dict['fasta_file_list_name'][1], 'w') as fasta_list_file:
        for fasta in fasta_list:
            fasta_list_file.write(f'{fasta}\n')
