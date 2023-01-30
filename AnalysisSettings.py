import numpy as np
import AccessiontoAlignment
import Analysis
from AccessiontoAlignment import alignment_finder
from sys import argv
from json import load
from os import listdir
from pathlib import Path

protein_info=argv[1]
argument_json=str(argv[2])
sequence_of_interest_fasta=argv[3]

with open(sequence_of_interest_fasta, 'r') as fasta:
    sequence_of_interest=''.join([x for x in fasta if x[0] != '>' if x != '']).strip().replace('\n', '')
with open(protein_info, 'r') as info_list:
    info_list=info_list.readlines()
    protein_list = [x.split()[-1] for x in info_list]
with open(argument_json, 'rb') as jfile:
    argument_dict=load(jfile)["arguments"]

character_to_replace=argument_dict['character_to_replace']
sequence_of_interest=[sequence_of_interest for x in protein_list]
msa_file=[argument_dict['msa_file_name'] for protein in protein_list]
native_plddts = [argument_dict['native_plddt'].replace(character_to_replace,protein) for protein in protein_list]
chimera_plddt=[argument_dict['chimera_plddt'].replace(character_to_replace,protein) for protein in protein_list]
plddt_files= native_plddts+chimera_plddt
number_of_subunits=[argument_dict['number_of_subunits'] for x in plddt_files]
reference_protein_name=[argument_dict['reference_protein_name'] for x in protein_list]
emboss_files=[argument_dict['emboss_names'].replace(character_to_replace,protein) for protein in protein_list]

if argument_dict['alphafold_output_folder'][0]!='#':
    alphafold_folders=[argument_dict['alphafold_output_folder'][1]+Path(file).stem for file in plddt_files]
    list(map(Analysis.generate_alphafold_files,alphafold_folders,plddt_files))
    reference_plddt=Analysis.generate_alphafold_files(argument_dict['alphafold_output_folder'][1]+argument_dict['reference_alphafold_folder'],
                                      argument_dict['reference_plddt'])

if number_of_subunits!=1 and argument_dict['averaged_native_plddt']!='#':
    averaged_native_plddt = [argument_dict['averaged_native_plddt'][1].replace(character_to_replace,protein) for protein in protein_list]
    averaged_chimera_plddt = [argument_dict['averaged_chimera_plddt'].replace(character_to_replace,protein) for protein in protein_list]
    averaged_plddt_files = native_plddts + chimera_plddt
    list(map(Analysis.averaging_multimer_plddt,plddt_files,averaged_plddt_files,number_of_subunits))
    averaged_reference_plddt=Analysis.averaging_multimer_plddt(argument_dict['reference_plddt'],argument_dict['averaged_reference'],number_of_subunits[0])
    plddt_files=averaged_plddt_files

if argument_dict['emboss_command'][0]!='#':
    splice_info = list(map(AccessiontoAlignment.alignment_finder, msa_file, sequence_of_interest, protein_list,
                           reference_protein_name, [argument_dict['emboss_command'][1] for x in emboss_files], emboss_files))
    splice_info=[x[1] for x in splice_info]
else:
    splice_info = list(map(AccessiontoAlignment.alignment_finder, msa_file, sequence_of_interest, protein_list,
                           reference_protein_name))
    splice_info=[x[1] for x in splice_info]

sequence_similarity=list(map(Analysis.get_sequence_similarity, emboss_files))

reference_boundaries=Analysis.get_reference_boundaries(sequence_of_interest[0],msa_file[0],argument_dict['reference_protein_name'])
if number_of_subunits==1:
    average_relative_stability=list(map(Analysis.relative_stability_full_protein,native_plddts,splice_info,chimera_plddt,
                                        reference_plddt,reference_boundaries))
else:
    average_relative_stability=list(map(Analysis.relative_stability_full_protein,averaged_native_plddt,splice_info,averaged_chimera_plddt,
                                        averaged_reference_plddt,reference_boundaries))

OverallDiff=list(map(Analysis.overall_confidence, native_plddts))
OverallChiDiff = list(map(Analysis.overall_confidence, chimera_plddt))

data_array=np.empty((len(protein_list) + 1, 5), dtype=object)
data_array[0,0], data_array[1:, 0]= 'Protein', protein_list
data_array[0,1], data_array[1:, 1]= 'S1 Sequence Similarity (%)',sequence_similarity
data_array[0,2], data_array[1:, 2]= 'Overall native plddt',OverallDiff
data_array[0,3], data_array[1:, 3]= 'Overall chimera plddt',OverallChiDiff
data_array[0,4], data_array[1:, 4]= 'Average Stability Difference',average_relative_stability

np.savetxt(argument_dict['analysis_output'], data_array, fmt="%s,%s,%s,%s,%s", delimiter="")
