from sys import argv
from json import load
from pathlib import Path
from numpy import savetxt, empty, delete
import Analysis
from AccessiontoAlignment import alignment_finder

argument_json = argv[1]
protein_info = argv[2]

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
character_to_replace = argument_dict['character_to_replace']
constant_sequence_of_interest=[constant_sequence_of_interest for name in comparison_proteins]
variant_sequence_of_interest=[variant_sequence_of_interest for name in comparison_proteins]
msa_file = [argument_dict['msa_file_name'] for protein in comparison_proteins]
native_plddts = [argument_dict['native_plddt'][1].replace(character_to_replace,protein) for protein in comparison_proteins]
chimera_plddt = [argument_dict['chimera_plddt'].replace(character_to_replace,protein) for protein in comparison_proteins]
plddt_files = native_plddts+chimera_plddt
number_of_subunits=[argument_dict['number_of_subunits'] for x in plddt_files]
constant_protein_name=[argument_dict['constant_fasta_identifier'] for x in comparison_proteins]
variant_protein_name=[argument_dict['variant_fasta_identifier'] for x in comparison_proteins]
constant_plddt = [argument_dict['constant_plddt'] for x in chimera_plddt]
emboss_files=[argument_dict['emboss_names'][1].replace(character_to_replace,protein) for protein in comparison_proteins]

if argument_dict['native_plddt'][0]=='':
    alphafold_folders=[argument_dict['alphafold_outputs_directory']+Path(file).stem+'/' for file in plddt_files]
    list(map(Analysis.generate_alphafold_files,alphafold_folders,plddt_files))
    Analysis.generate_alphafold_files(argument_dict['alphafold_outputs_directory'] + argument_dict['constant_alphafold_folder_name'] + '/',
                                      argument_dict['constant_plddt'])
if number_of_subunits!=1:
    averaged_native_plddt = [argument_dict['averaged_native_plddt'].replace(character_to_replace,protein) for protein in comparison_proteins]
    averaged_chimera_plddt = [argument_dict['averaged_chimera_plddt'].replace(character_to_replace,protein) for protein in comparison_proteins]
    list(map(Analysis.averaging_multimer_plddt,native_plddts,averaged_native_plddt,number_of_subunits))
    list(map(Analysis.averaging_multimer_plddt, chimera_plddt, averaged_chimera_plddt, number_of_subunits))
    Analysis.averaging_multimer_plddt(argument_dict['constant_plddt'], argument_dict['averaged_constant'],
                                      number_of_subunits[0])
    averaged_constant_plddt=[argument_dict['averaged_constant'] for x in averaged_chimera_plddt]
constant_boundary=alignment_finder(constant_fasta, constant_sequence_of_interest[0], constant_protein_name[0], constant_protein_name[0])[1]

if argument_dict['emboss_command'][0]!='' or argument_dict['emboss_names'][0]!='':
    variant_splice_info = list(map(alignment_finder, msa_file, variant_sequence_of_interest, comparison_proteins,
                                   variant_protein_name))
    variant_splice_info = [x[1] for x in variant_splice_info]
else:
    variant_splice_info = list(map(alignment_finder, msa_file, variant_sequence_of_interest, comparison_proteins,
                                   variant_protein_name, [argument_dict['emboss_command'][1] for x in emboss_files], emboss_files))
    variant_splice_info=[x[1] for x in variant_splice_info]
# MAKE SURE YOU KNOW WHICH PROTEIN GOES FIRST
average_relative_stability=[]
if number_of_subunits==1:
    for index,boundary in enumerate(variant_splice_info):
        if argument_dict['constant_protein'][0]=='#' and argument_dict['variant_protein'][0]=='':
            chimera_boundary=[(0,len(constant_sequence_of_interest[0])),(len(constant_sequence_of_interest[0]),None)]
            constant_stability = Analysis.relative_stability(constant_plddt[index],constant_boundary,chimera_plddt[index],chimera_boundary[0])
            variant_stability=Analysis.relative_stability(native_plddts[index],variant_splice_info[index],chimera_plddt[index],chimera_boundary[1])
            average_relative_stability.append((constant_stability[0]+variant_stability[0])/(constant_stability[1]+variant_stability[1]))

        if argument_dict['variant_protein'][0]=='#' and argument_dict['constant_protein'][0]=='':
            chimera_boundary = [(0, len(constant_sequence_of_interest[0])),
                                (len(constant_sequence_of_interest[0]), None)]
            variant_stability = Analysis.relative_stability(native_plddts[index], variant_splice_info[index],
                                                            chimera_plddt[index], chimera_boundary[0])
            constant_stability = Analysis.relative_stability(constant_plddt[index], constant_boundary, chimera_plddt[index],
                                                             chimera_boundary[1])
            average_relative_stability.append(
                (constant_stability[0] + variant_stability[0]) / (constant_stability[1] + variant_stability[1]))
        
else:
    for index, boundary in enumerate(variant_splice_info):
        if argument_dict['constant_protein'][0] == '#' and argument_dict['variant_protein'][0] == '':
            chimera_boundary = [(0, len(constant_sequence_of_interest[0])),
                                (len(constant_sequence_of_interest[0]), None)]
            constant_stability = Analysis.relative_stability(averaged_constant_plddt[index], constant_boundary, averaged_chimera_plddt[index],
                                                             chimera_boundary[0])
            variant_stability = Analysis.relative_stability(averaged_native_plddt[index], variant_splice_info[index],
                                                            averaged_chimera_plddt[index], chimera_boundary[1])
            average_relative_stability.append(
                (constant_stability[0] + variant_stability[0]) / (constant_stability[1] + variant_stability[1]))
            print(averaged_constant_plddt[index], constant_boundary, averaged_chimera_plddt[index],
                  chimera_boundary[0])
            print(averaged_native_plddt[index], variant_splice_info[index],
                  averaged_chimera_plddt[index], chimera_boundary[1])
        if argument_dict['variant_protein'][0] == '#' and argument_dict['constant_protein'][0] == '':
            chimera_boundary = [(0, len(constant_sequence_of_interest[0])),
                                (len(constant_sequence_of_interest[0]), None)]
            variant_stability = Analysis.relative_stability(averaged_native_plddt[index], variant_splice_info[index],
                                                            averaged_chimera_plddt[index], chimera_boundary[0])
            constant_stability = Analysis.relative_stability(averaged_constant_plddt[index], constant_boundary, averaged_chimera_plddt[index],
                                                             chimera_boundary[1])
            average_relative_stability.append(
                (constant_stability[0] + variant_stability[0]) / (constant_stability[1] + variant_stability[1]))

# Absolute stability which averages plddt scores across all residue positions will be calculated for the native splice partners
# and their chimeras
overall_stability = list(map(Analysis.overall_confidence, native_plddts))
overall_chimera_stability = list(map(Analysis.overall_confidence, chimera_plddt))
# Column names for each piece of data calculated are recorded from the json, column names with a # in front as described at
# the top will not be included in the array
# DATA RECORDED OUTSIDE OF THESE FIVE OBJECTS IS NOT CUSTOMIZABLE NO MATTER THE JSON FILE ARGUMENTS, THEY CAN BE TURNED ON AND OFF
# BUT NOT SWITCHED WITHOUT PERSONAL EDITING OF THIS CODE
data_array = empty((len(comparison_proteins) + 1, 5), dtype=object)
column_names = [x[1] for x in argument_dict['analysis_column_names']]
preferred_columns = [x[0] for x in argument_dict['analysis_column_names']]
# If a # is indicated in front of the emboss_names it will be excluded from the data output file
if argument_dict['emboss_names'][0] == '':
    sequence_similarity=list(map(Analysis.get_sequence_similarity, emboss_files))
else:
    sequence_similarity = ['' for x in comparison_proteins]
    preferred_columns[1] = '#'
types_of_data = [comparison_proteins, sequence_similarity, overall_stability, overall_chimera_stability, average_relative_stability]
column_count = 0
# Column names with # in front are excluded from the final array, otherwise the name and accompanying data are recorded into the array
for name, corresponding_data, preferred in zip(column_names, types_of_data, preferred_columns):
    if preferred=='':
        data_array[0, column_count], data_array[1:, column_count] = name, corresponding_data
        column_count+=1
    else:
        data_array=delete(data_array, -1, 1)
savetxt(argument_dict['analysis_output_csv'], data_array, fmt=','.join('%s' for x in preferred_columns if x==''), delimiter=",")
