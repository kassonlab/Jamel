def generate_alphafold_files(output_folder, new_plddt='NA',new_pdb='NA'):
    """Creates a text file containing the plddt values of the highest ranked model, extracted from alphafold's result pkl file."""
    from pickle import load as pload
    from  json import load as jload
    from os import path
    from shutil import copy
    from numpy import savetxt
    # Checking to see if ranking_debug.json exists. This file is the last to be output by alphafold and is a check that
    # the pkl file you want to extract from exists, as well as to avoid errors
    if path.exists(output_folder + 'ranking_debug.json'):
        if new_pdb!= 'NA':
            copy(output_folder + 'ranked_0.pdb', new_pdb)
        if new_plddt!= 'NA':
            # ranking_debug is also useful for determining which result pkl file is highest ranked. The model result pkls are
            # numbered by the order they are created and not their overall confidence score. The information about their rank by
            # score is found in ranking_debug.json
            with open(output_folder + 'ranking_debug.json', 'r') as jfile:
                HighestRankModel=jload(jfile)['order'][0]
                with open({HighestRankModel}+'.pkl', 'rb') as pfile:
                    data=pload(pfile)
                    savetxt(new_plddt, data['plddt'], fmt="%s", delimiter=" ")


def get_sequence_similarity(emboss_file):
    """Returns sequence similarity from an emboss needle file."""
    with open(emboss_file,'r') as file:
        file=file.read().split('#')
        for line in file:
            if 'Similarity' in line:
                emboss_score=line.split()[-1].replace('(','').replace(')','').replace('%','')
    return emboss_score


def overall_confidence(plddt_file):
    """Returns the average confidence score from a protein's plddt file."""
    plddt= [float(score) for score in open(plddt_file, 'r').readlines()]
    average_plddt=sum(plddt)/len(plddt)
    return average_plddt


def get_reference_boundaries(sequence_of_interest, msa, fasta_identifier):
    with open(msa, "r") as alignment:
        alignment = alignment.read().split('>')
        sequence_dictionary = {sequence.split('\n')[0]: ''.join(sequence.split('\n')[1:]) for sequence in alignment if
                               len(sequence) != 0}
    reference_sequence = ''.join([x for x in sequence_dictionary[fasta_identifier] if x != '-'])

    splice_start = reference_sequence.find(sequence_of_interest)
    splice_end = splice_start + len(sequence_of_interest)
    boundaries = list({0, splice_start, splice_end, len(reference_sequence) - 1})
    boundaries.sort()
    spliced_out = (splice_start, splice_end)
    boundaries_tuples = []
    for x in range(len(boundaries) - 1):
        if (boundaries[x], boundaries[x + 1]) == spliced_out:
            boundaries_tuples.append(('NS', boundaries[x], boundaries[x + 1]))
        else:
            boundaries_tuples.append(('S', boundaries[x], boundaries[x + 1]))
    return boundaries_tuples


def relative_stability(native_plddt, native_boundary_tuple,chimera_plddt, chimera_boundary_tuple):
    """Returns the relative percent difference between the two equally sized sections of plddt scores that are outlined with
    native_boundary_tuple and chimera_boundary_tuple. relative_difference=(compared value-reference value)/reference value * 100
    Native scores are assumed to be the reference value in this formula for relative difference"""
    # Pulling the plddt values as floats that start at native_boundary_tuple[0] and chimera_boundary_tuple[0], and end at
    # native_boundary_tuple[1] and chimera_boundary_tuple[1] but dont include index [1] scores.
    native_score = [float(score) for score in open(native_plddt, 'r').readlines()][native_boundary_tuple[0]:native_boundary_tuple[1]]
    chimera_score=[float(score) for score in open(chimera_plddt, 'r').readlines()][chimera_boundary_tuple[0]:chimera_boundary_tuple[1]]
    # Recording the length of the residue scores for averaging purposes later
    splice_length=len(chimera_score)
    relative_difference=sum([(chimera-native)/native*100 for native,chimera in zip(native_score,chimera_score)])
    return relative_difference,splice_length


def relative_stability_full_protein(native_plddt, native_boundary_tuple,chimera_plddt,reference_plddt,reference_tuples):
    native_score = [float(score) for score in open(native_plddt, 'r').readlines()]
    chimera_score = [float(score) for score in open(chimera_plddt, 'r').readlines()]
    raw_stability = 0
    current_chimera_index=0
    for index,tuples in enumerate(reference_tuples):
        if tuples[0] == 'NS':
            comparison_splice_length = native_boundary_tuple[1] - native_boundary_tuple[0]
            raw_stability += relative_stability(native_score, native_boundary_tuple, chimera_score,
                                                (current_chimera_index, current_chimera_index + comparison_splice_length))[0]
            current_chimera_index+=comparison_splice_length
        elif tuples[0] == 'S':
            reference_splice_length=tuples[1]-tuples[0]
            raw_stability += relative_stability(reference_plddt,tuples[1:],chimera_score,(current_chimera_index,current_chimera_index+reference_splice_length))[0]
            current_chimera_index += reference_splice_length
    averaged_relative_stability=raw_stability/len(chimera_plddt)
    return averaged_relative_stability


def averaging_multimer_plddt(plddt_file, new_plddt_file,subunits):
    """This function takes a plddt and averages the scores
    for each residue position across the number of subunints specified"""
    # Using list comprehension to turn the plddt file into a list of floats
    multimer_plddt=[float(score) for score in open(plddt_file, 'r').readlines()]
    # Calculating the length a subunits to have for step size when iterating through the list later
    monomer_length=int(len(multimer_plddt) / int(subunits))
    # creating a file to input the averaged scores
    new_plddt = open(new_plddt_file, 'w')
    # using list comprehension to step through each the residue position of each subunit and
    # collect their scores, average them and return them to the new list
    averaged_scores=[sum(multimer_plddt[residue_index::monomer_length])/subunits for residue_index in range(monomer_length)]
    # Looping through the new list and inputting the averaged scores into the new file that was created
    for score in averaged_scores:
        new_plddt.write(f'{score}\n')
    new_plddt.close()




