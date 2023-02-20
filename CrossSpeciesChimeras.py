import ChimeraGenerator
from sys import argv
from json import load
import AccessiontoAlignment
argument_json=argv[1]
protein_info=argv[2]


with open(argument_json, 'rb') as jfile:
    argument_dict = load(jfile)["arguments"]
with open(protein_info, 'r') as info_list:
    info_list = info_list.readlines()
    accession_number = [x.split()[0] for x in info_list]
    comparison_proteins = [x.split()[1] for x in info_list]
constant_sequence_of_interest=argument_dict['constant_sequence_of_interest']
variant_sequence_of_interest=argument_dict['variant_sequence_of_interest']
with open(constant_sequence_of_interest, 'r') as fasta:
    constant_sequence_of_interest = ''.join([x for x in fasta if x[0] != '>' if x != '']).strip().replace('\n', '')
with open(variant_sequence_of_interest, 'r') as fasta:
    variant_sequence_of_interest = ''.join([x for x in fasta if x[0] != '>' if x != '']).strip().replace('\n', '')

constant_fasta=argument_dict['constant_protein'][1]
variant_fasta=argument_dict['variant_protein'][1]
character_to_replace=argument_dict['character_to_replace']
constant_fasta_identifier=argument_dict['constant_fasta_identifier']
variant_fasta_identifier=[argument_dict['variant_fasta_identifier'] for name in comparison_proteins]
email=[argument_dict['email_for_accession'] for name in comparison_proteins]
subunits=[argument_dict['number_of_subunits'] for name in comparison_proteins]
msa=[argument_dict['msa_file_name'] for name in comparison_proteins]
constant_sequence_of_interest=[constant_sequence_of_interest for name in comparison_proteins]
variant_sequence_of_interest=[variant_sequence_of_interest for name in comparison_proteins]
monomer_fastas=[argument_dict['monomer_fastas'][1].replace(character_to_replace,name) for name in comparison_proteins]
multimer_fastas=[argument_dict['multimer_fastas'][1].replace(character_to_replace,name) for name in comparison_proteins]
chimera_fastas=[argument_dict['chimera_fastas'].replace(character_to_replace,name) for name in comparison_proteins]

if subunits[0]==1 and argument_dict['monomer_fastas'][0]=='':
    list(map(AccessiontoAlignment.accession_to_fasta, monomer_fastas, accession_number, email, subunits))
elif subunits[0]>1 and argument_dict['multimer_fastas'][0]=='':
    list(map(AccessiontoAlignment.accession_to_fasta,monomer_fastas, accession_number, email,subunits, multimer_fastas))

if argument_dict['muscle_command_for_msa'][0]=='':
    AccessiontoAlignment.multiple_sequence_alignment(monomer_fastas, argument_dict['msa_fasta'], msa[0],
                                                     variant_fasta, argument_dict['muscle_command_for_msa'][1])
variant_sequences=[x[0] for x in map(AccessiontoAlignment.alignment_finder, msa, variant_sequence_of_interest, comparison_proteins,
                                     variant_fasta_identifier)]
if argument_dict['constant_protein'][0]=='#' and argument_dict['variant_protein'][0]=='':
    chimera_sequences=[constant_sequence_of_interest[0] + variants for variants in variant_sequences]
elif argument_dict['constant_protein'][0]=='' and argument_dict['variant_protein'][0]=='#':
    chimera_sequences=[variants + constant_sequence_of_interest[0] for variants in variant_sequences]
list(map(ChimeraGenerator.fasta_creation,chimera_fastas,chimera_sequences,subunits))
if argument_dict['fasta_file_list_name'][0] == '':
    if subunits == 1:
        fasta_list = monomer_fastas + chimera_fastas
    else:
        fasta_list = multimer_fastas + chimera_fastas
    with open(argument_dict['fasta_file_list_name'][1], 'w') as fasta_list_file:
        for fasta in fasta_list:
            fasta_list_file.write(f'{fasta}\n')
