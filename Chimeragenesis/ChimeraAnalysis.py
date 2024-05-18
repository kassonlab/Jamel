#!python
#
# Code 2023 by Jamel Simpson

""" This script is made for the analysis of chimeric proteins between one designated protein called a reference protein,
and a designated section within the reference protein sequence, called sequence_of_interest being replaced by homologous sequences outlined in a list
accession number sequences in a protein_info file"""

from sys import argv
from json import load
from pathlib import Path
from numpy import savetxt, empty, delete
import Analysis
from AccessiontoAlignment import alignment_finder



# 3 command line inputs are required for the analysis script, similarly to the production.script, a .tsv file with accession numbers in column one
# and the corresponding name you want associated with that sequence, these name will be used throughout all naming
# including:fasta file names, alphafold out folders, in file fasta identifiers, plddt data files and so on
# DO NOT INCLUDE YOUR REFERENCE PROTEIN IN THIS FILE
protein_info = argv[1]
# The second is a analysis_arguments.json file modified to your liking, this file has naming conventions for all important
# files generated like plddt files, emboss files and the output file containing the tables of data associated with your
# chimeras, keep naming conventions consistent from chimera_arguments.json for the basename of outputs, the naming works by replacing the
# character_to_replace which is defaulted with an asterisk ('*'), with the names supplied in the protein_info list that was generated
# for the last input, there are also a couple operation that come as arrays with enclosed in [] and the first input is just
# quotations "", if you place a # within the quotes it will prevent that file from being created, for example you prevent
# the emboss files from being included in the data analysis by the script using this method like this ["#","/gpfs/gpfs0/scratch/jws6pq/Notebook/Emboss/*S1vsSARS2.emboss"] or
# just include already created emboss files but not run emboss["#","/scratch/jws6pq/EMBOSS-6.6.0/emboss/needle"] or
# preventing certain data columns from your final output like [["","Protein"],["","Sequence Similarity (%)"],["#","Overall Native Score"],["","Overall Chimera Score"],["","Relative Stability"]]
# which would delete Overall Native Score and data
# PROVIDE YOUR REFERENCE IN THE JSON INPUT FILE
argument_json = str(argv[2])
# Lastly and most simply, create a fasta file with the sequence that you want spliced out and replaced of your reference protein
sequence_of_interest_fasta = argv[3]
# These lines are extracting all the inputs and info outlined in the previous comments
with open(sequence_of_interest_fasta, 'r') as fasta:
    sequence_of_interest = ''.join([x for x in fasta if x[0] != '>' if x != '']).strip().replace('\n', '')
with open(protein_info, 'r') as info_list:
    info_list = info_list.readlines()
    protein_list = [x.split()[-1] for x in info_list]
with open(argument_json, 'rb') as jfile:
    argument_dict = load(jfile)["arguments"]
# Here, important information from your command line inputs are turned into list iterables to be used in the map function later
character_to_replace = argument_dict['character_to_replace']
sequence_of_interest = [sequence_of_interest for x in protein_list]
msa_file = [argument_dict['msa_file_name'] for protein in protein_list]
native_plddts = [argument_dict['native_plddt'][1].replace(character_to_replace,protein) for protein in protein_list]
chimera_plddt = [argument_dict['chimera_plddt'].replace(character_to_replace,protein) for protein in protein_list]
plddt_files = native_plddts+chimera_plddt
number_of_subunits=[argument_dict['number_of_subunits'] for x in plddt_files]
reference_protein_name=[argument_dict['reference_protein_name'] for x in protein_list]
emboss_files=[argument_dict['emboss_names'][1].replace(character_to_replace,protein) for protein in protein_list]
# If no # was indicated next to alphafold_output_folder argument in json, the script will create plddt files for all proteins
# in protein_info including the reference and chimeras
# CONSISTENT NAMING IS IMPORTANT, ALPHAFOLD FOLDERS ARE NAMED AFTER FASTA FILES, MAKE SURE PLDDT AND FASTA HAVE EXACT
# SAME NAMING CONVENTIONS AND THAT ALL NECESSARY ALPHAFOLD FOLDERS ARE IN THE SAME PLACE
if argument_dict['native_plddt'][0]=='':
    alphafold_folders=[argument_dict['alphafold_outputs_directory']+Path(file).stem+'/' for file in plddt_files]
    list(map(Analysis.generate_alphafold_files,alphafold_folders,plddt_files))
    Analysis.generate_alphafold_files(argument_dict['alphafold_outputs_directory'] + argument_dict['reference_alphafold_folder_name'] + '/',
                                      argument_dict['reference_plddt'])
    reference_plddt=[argument_dict['reference_plddt'] for x in chimera_plddt]
# If no # was indicated in front of the averaged_native_plddt argument and you have a multimer protein, all plddt files will be averaged across
# subunits to a single score for every residue position
if number_of_subunits!=1 and argument_dict['averaged_native_plddt'][0]=="":
    averaged_native_plddt = [argument_dict['averaged_native_plddt'][1].replace(character_to_replace,protein) for protein in protein_list]
    averaged_chimera_plddt = [argument_dict['averaged_chimera_plddt'].replace(character_to_replace,protein) for protein in protein_list]
    list(map(Analysis.averaging_multimer_plddt,native_plddts,averaged_native_plddt,number_of_subunits))
    list(map(Analysis.averaging_multimer_plddt, chimera_plddt, averaged_chimera_plddt, number_of_subunits))
    Analysis.averaging_multimer_plddt(argument_dict['reference_plddt'], argument_dict['averaged_reference'],
                                      number_of_subunits[0])
    averaged_reference_plddt=[argument_dict['averaged_reference'] for x in averaged_chimera_plddt]
# If no # was indicated before the emboss_command or the emboss_names, emboss files calculating the sequence similarity between
# the sequence_of_interest and homologous sequences from protein outlined in protein_info will be created
# No matter what, information will be put into splice_info that indicates the boundaries for homologous sequences
# previously described within the native splice partners
if argument_dict['emboss_command'][0]!='' or argument_dict['emboss_names'][0]!='':
    splice_info = list(map(alignment_finder, msa_file, sequence_of_interest, protein_list,
                           reference_protein_name))
    splice_info = [x[1] for x in splice_info]
else:
    splice_info = list(map(alignment_finder, msa_file, sequence_of_interest, protein_list,
                           reference_protein_name, [argument_dict['emboss_command'][1] for x in emboss_files], emboss_files))
    splice_info=[x[1] for x in splice_info]
# Averaged relative stability for all chimeras are calculated
if number_of_subunits==1:
    average_relative_stability=list(map(Analysis.average_relative_stability_full_chimera,native_plddts,splice_info,chimera_plddt,
                                        reference_plddt,sequence_of_interest,msa_file,reference_protein_name))
else:
    average_relative_stability=list(map(Analysis.average_relative_stability_full_chimera,averaged_native_plddt,splice_info,averaged_chimera_plddt,
                                        averaged_reference_plddt,sequence_of_interest,msa_file,reference_protein_name))
# Absolute stability which averages plddt scores across all residue positions will be calculated for the native splice partners
# and their chimeras
overall_stability = list(map(Analysis.overall_confidence, native_plddts))
overall_chimera_stability = list(map(Analysis.overall_confidence, chimera_plddt))
# Column names for each piece of data calculated are recorded from the json, column names with a # in front as described at
# the top will not be included in the array
# DATA RECORDED OUTSIDE OF THESE FIVE OBJECTS IS NOT CUSTOMIZABLE NO MATTER THE JSON FILE ARGUMENTS, THEY CAN BE TURNED ON AND OFF
# BUT NOT SWITCHED WITHOUT PERSONAL EDITING OF THIS CODE
data_array = empty((len(protein_list) + 1, 5), dtype=object)
column_names = [x[1] for x in argument_dict['analysis_column_names']]
preferred_columns = [x[0] for x in argument_dict['analysis_column_names']]
# If a # is indicated in front of the emboss_names it will be excluded from the data output file
if argument_dict['emboss_names'][0] == '':
    sequence_similarity=list(map(Analysis.get_sequence_similarity, emboss_files))
else:
    sequence_similarity = ['' for x in protein_list]
    preferred_columns[1] = '#'
types_of_data = [protein_list, sequence_similarity, overall_stability, overall_chimera_stability, average_relative_stability]
column_count = 0
# Column names with # in front are excluded from the final array, otherwise the name and accompanying data are recorded into the array
for name, corresponding_data, preferred in zip(column_names, types_of_data, preferred_columns):
    if preferred=='':
        data_array[0, column_count], data_array[1:, column_count] = name, corresponding_data
        column_count+=1
    else:
        data_array=delete(data_array, -1, 1)
savetxt(argument_dict['analysis_output_csv'], data_array, fmt=','.join('%s' for x in preferred_columns if x==''), delimiter=",")
