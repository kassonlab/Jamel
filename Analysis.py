#! python
#
# Code 2023 by Jamel Simpson


from pickle import load as p_load
from json import load as j_load
from os import system, path
from shutil import copy
from numpy import savetxt

def run_Foldx(foldx_file,pdb_file,foldx_command):
    """Call FoldX.  This is optional functionality."""
    pdb_dir = path.dirname(pdb_file)
    foldx_dir = path.dirname(foldx_file)
    system(f'{foldx_command}  -c Stability --pdb {path.basename(pdb_file)} '
           f'--output-dir {foldx_dir} --output-file {path.basename(foldx_file)} '
           f'--pdb-dir {pdb_dir}')

def get_Foldx_results(foldx_file):
    """Read FoldX results."""
    with open(foldx_file, 'r') as foldx_score:
        return foldx_score.read().split()[1]

def generate_alphafold_files(output_folder, new_plddt='NA', new_pdb='NA'):
    """Creates a text file containing the plddt values of the highest_rank_model extracted from alphafold's result pkl file
    and renames the ranked_0.pdb file and places it in the desired directory."""
    # Checking to see if ranking_debug.json exists. This file is the last to be output by alphafold and is a check that
    # the pkl file you want to extract from exists, as well as to avoid errors
    if path.exists(output_folder + 'ranking_debug.json'):
        if new_pdb != 'NA':
            # The highest ranked structure is copied with a new name and directory
            copy(output_folder + 'ranked_0.pdb', new_pdb)
        if new_plddt != 'NA':
            # ranking_debug is also useful for determining which result pkl file is the highest ranked. The model result pkl files are
            # numbered by the order they are created and not their overall confidence score. The information about their rank by
            # confidence score is found in ranking_debug.json
            with open(output_folder + 'ranking_debug.json', 'r') as jfile:
                highest_rank_model = j_load(jfile)['order'][0]
                with open(f'{output_folder}result_{highest_rank_model}.pkl', 'rb') as pfile:
                    data = p_load(pfile)
                    # The plddt scores are put into a column in a text file named by new_plddt
                    savetxt(new_plddt, data['plddt'], fmt='%s', delimiter=' ')

def get_sequence_similarity(emboss_file):
    """Returns sequence similarity from an emboss needle file."""
    with open(emboss_file, 'r') as infile:
        infile = infile.read().split('#')
        for line in infile:
            if 'Similarity' in line:
                emboss_score = line.split()[-1].replace('(','').replace(')','').replace('%','')
    return emboss_score

def overall_confidence(plddt_file):
    """Returns the average confidence score from a protein's plddt file."""
    with open(plddt_file, 'r') as infile:
        plddt = tuple(float(score) for score in infile.readlines())
    average_plddt = sum(plddt)/len(plddt)
    return average_plddt

def get_reference_boundaries(sequence_of_interest, msa, fasta_identifier):
    """Returns the list_of_boundary_tuples within the reference protein that contain the sequence_of_interest, as well as, the boundaries of the
    sections before and after. Boundary tuples that contain the sequence_of_interest are marked by 'NS' as in Not Spliced into
    the resulting chimera, and all other tuples are marked with 'S' as in spliced into the chimera.
    The tuples are provided as so: ('NS',boundary_one,boundary_two) or ('S',boundary_one,boundary_two)"""
    # This uses the msa to grab the reference sequence outlined by fasta_identifier
    with open(msa, 'r') as alignment:
        alignment = alignment.read().split('>')
        sequence_dictionary = {sequence.split('\n')[0]: ''.join(sequence.split('\n')[1:]) for sequence in alignment if
                               len(sequence) != 0}
    reference_sequence = ''.join(x for x in sequence_dictionary[fasta_identifier] if x != '-')
    # This is recording the 'NS' boundaries that indicate the boundaries of the sequence_of_interest
    splice_start = reference_sequence.find(sequence_of_interest)
    splice_end = splice_start + len(sequence_of_interest)
    # Then those boundaries are compared against the very beginning and end of the proteins, by introducing them into a set
    # to get of redundancy if the sequence_of_interest boundaries contain the beginning or end
    boundaries = tuple({0, splice_start, splice_end, len(reference_sequence)})
    # They are sorted into ascending order
    boundaries=sorted(boundaries)
    spliced_out = (splice_start, splice_end)
    list_of_boundary_tuples = []
    # Then the loop checks if they're the sequence_of_interest boundaries and marks them accordingly
    for x in range(len(boundaries) - 1):
        if (boundaries[x], boundaries[x + 1]) == spliced_out:
            list_of_boundary_tuples.append(('NS', boundaries[x], boundaries[x + 1]))
        else:
            list_of_boundary_tuples.append(('S', boundaries[x], boundaries[x + 1]))
    return list_of_boundary_tuples

def relative_stability(native_plddt, native_boundary_tuple, chimera_plddt, chimera_boundary_tuple):
    """Returns the relative percent difference between the two equally sized sections of plddt scores that are outlined with
    native_boundary_tuple and chimera_boundary_tuple. relative_difference=(compared value-reference value)/reference value * 100
    Native scores are assumed to be the reference value in this formula for relative difference"""
    # Pulling the plddt values as floats that start at native_boundary_tuple[0] and chimera_boundary_tuple[0], and end at
    # native_boundary_tuple[1] and chimera_boundary_tuple[1] but dont include index [1] scores.
    with open(native_plddt, 'r') as infile:
        native_score = tuple(float(score) for
                        score in infile.readlines())[native_boundary_tuple[0]:native_boundary_tuple[1]]
    with open(chimera_plddt, 'r') as infile:
        chimera_score = tuple(float(score) for
                         score in infile.readlines())[chimera_boundary_tuple[0]:chimera_boundary_tuple[1]]
    # Recording the length of the residue scores for averaging purposes later
    splice_length = len(chimera_score)
    relative_difference = sum((chimera-native) / native*100 for native, chimera in zip(native_score, chimera_score))
    return relative_difference, splice_length

def average_relative_stability_full_chimera(native_plddt, native_boundary_tuple,
                                            chimera_plddt, reference_plddt,
                                            sequence_of_interest, msa, reference_fasta_identifier):
    """Returns the averaged relative stability of a full length chimera, given a: multiple sequence alignment, reference
    protein's plddt and identifier in the msa, the plddt of the wild-type splice partner for the reference and the boundaries as a tuple
    that contain the spliced in sequence, and the resulting chimera's plddt"""
    raw_stability = 0
    # This variable is recording the chimera boundaries for relative stability comparison and will also help withe averaging later
    current_chimera_index = 0
    # Retrieves boundaries for which sections of the reference protein to compare to the chimera
    reference_boundaries = get_reference_boundaries(sequence_of_interest,msa,reference_fasta_identifier)
    for index, tuples in enumerate(reference_boundaries):
        # NS indicates that its time to calculate relative stability against the parent splice partner rather than the
        #reference protein
        if tuples[0] == 'NS':
            comparison_splice_length = native_boundary_tuple[1] - native_boundary_tuple[0]
            raw_stability += relative_stability(native_plddt, native_boundary_tuple, chimera_plddt,
                                                (current_chimera_index, current_chimera_index + comparison_splice_length))[0]
            current_chimera_index += comparison_splice_length
        # S represents the opposite, that the chimera should now be compared to the reference protein
        elif tuples[0] == 'S':
            reference_splice_length = tuples[2] - tuples[1]
            raw_stability += relative_stability(reference_plddt,
                                                tuples[1:],
                                                chimera_plddt,
                                                (current_chimera_index,
                                                 current_chimera_index + reference_splice_length))[0]
            current_chimera_index += reference_splice_length
    averaged_relative_stability = raw_stability / (current_chimera_index)
    return averaged_relative_stability

def averaging_multimer_plddt(plddt_file, new_plddt_file,subunits):
    """This function takes a plddt and averages the scores
    for each residue position across the number of subunints specified"""
    # Using list comprehension to turn the plddt file into a list of floats
    with open(plddt_file, 'r') as infile:
        multimer_plddt = tuple(float(score) for score in infile.readlines())
    # Calculating the length a subunits to have for step size when iterating through the list later
    monomer_length = int(len(multimer_plddt) / int(subunits))
    # creating a file to input the averaged scores
    with open(new_plddt_file, 'w') as new_plddt:
        # using list comprehension to step through each the residue position of each subunit and
        # collect their scores, average them and return them to the new list
        averaged_scores = tuple(sum(multimer_plddt[residue_index::monomer_length])/subunits
                           for residue_index in range(monomer_length))
        # Looping through the new list and inputting the averaged scores into the new file that was created
        for score in averaged_scores:
            new_plddt.write(f'{score}\n')
