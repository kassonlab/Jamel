def blosum_62_matrix():
    """Simply returns the blosum 62 matrix"""
    from numpy import array
    blosum_62_array = array([['','A','R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W','Y', 'V', 'B', 'Z', 'X', '-'],
                           ['A', 4.0, -1.0, -2.0, -2.0, 0.0, -1.0, -1.0, 0.0, -2.0, -1.0, -1.0, -1.0, -1.0, -2.0, -1.0, 1.0,0.0, -3.0, -2.0, 0.0, -2.0, -1.0, 0.0, -4.0],
                           ['R', -1.0, 5.0, 0.0, -2.0, -3.0, 1.0, 0.0, -2.0, 0.0, -3.0, -2.0, 2.0, -1.0, -3.0, -2.0, -1.0,-1.0, -3.0, -2.0, -3.0, -1.0, 0.0, -1.0, -4.0],
                           ['N', -2.0, 0.0, 6.0, 1.0, -3.0, 0.0, 0.0, 0.0, 1.0, -3.0, -3.0, 0.0, -2.0, -3.0, -2.0, 1.0, 0.0,-4.0, -2.0, -3.0, 3.0, 0.0, -1.0, -4.0],
                           ['D', -2.0, -2.0, 1.0, 6.0, -3.0, 0.0, 2.0, -1.0, -1.0, -3.0, -4.0, -1.0, -3.0, -3.0, -1.0, 0.0,-1.0, -4.0, -3.0, -3.0, 4.0, 1.0, -1.0, -4.0],
                           ['C', 0.0, -3.0, -3.0, -3.0, 9.0, -3.0, -4.0, -3.0, -3.0, -1.0, -1.0, -3.0, -1.0, -2.0, -3.0,-1.0, -1.0, -2.0, -2.0, -1.0, -3.0, -3.0, -2.0, -4.0],
                           ['Q', -1.0, 1.0, 0.0, 0.0, -3.0, 5.0, 2.0, -2.0, 0.0, -3.0, -2.0, 1.0, 0.0, -3.0, -1.0, 0.0,-1.0, -2.0, -1.0, -2.0, 0.0, 3.0, -1.0, -4.0],
                           ['E', -1.0, 0.0, 0.0, 2.0, -4.0, 2.0, 5.0, -2.0, 0.0, -3.0, -3.0, 1.0, -2.0, -3.0, -1.0, 0.0,-1.0, -3.0, -2.0, -2.0, 1.0, 4.0, -1.0, -4.0],
                           ['G', 0.0, -2.0, 0.0, -1.0, -3.0, -2.0, -2.0, 6.0, -2.0, -4.0, -4.0, -2.0, -3.0, -3.0, -2.0, 0.0,-2.0, -2.0, -3.0, -3.0, -1.0, -2.0, -1.0, -4.0],
                           ['H', -2.0, 0.0, 1.0, -1.0, -3.0, 0.0, 0.0, -2.0, 8.0, -3.0, -3.0, -1.0, -2.0, -1.0, -2.0, -1.0,-2.0, -2.0, 2.0, -3.0, 0.0, 0.0, -1.0, -4.0],
                           ['I', -1.0, -3.0, -3.0, -3.0, -1.0, -3.0, -3.0, -4.0, -3.0, 4.0, 2.0, -3.0, 1.0, 0.0, -3.0, -2.0,-1.0, -3.0, -1.0, 3.0, -3.0, -3.0, -1.0, -4.0],
                           ['L', -1.0, -2.0, -3.0, -4.0, -1.0, -2.0, -3.0, -4.0, -3.0, 2.0, 4.0, -2.0, 2.0, 0.0, -3.0, -2.0,-1.0, -2.0, -1.0, 1.0, -4.0, -3.0, -1.0, -4.0],
                           ['K', -1.0, 2.0, 0.0, -1.0, -3.0, 1.0, 1.0, -2.0, -1.0, -3.0, -2.0, 5.0, -1.0, -3.0, -1.0, 0.0,-1.0, -3.0, -2.0, -2.0, 0.0, 1.0, -1.0, -4.0],
                           ['M', -1.0, -1.0, -2.0, -3.0, -1.0, 0.0, -2.0, -3.0, -2.0, 1.0, 2.0, -1.0, 5.0, 0.0, -2.0, -1.0,-1.0, -1.0, -1.0, 1.0, -3.0, -1.0, -1.0, -4.0],
                           ['F', -2.0, -3.0, -3.0, -3.0, -2.0, -3.0, -3.0, -3.0, -1.0, 0.0, 0.0, -3.0, 0.0, 6.0, -4.0, -2.0,-2.0, 1.0, 3.0, -1.0, -3.0, -3.0, -1.0, -4.0],
                           ['P', -1.0, -2.0, -2.0, -1.0, -3.0, -1.0, -1.0, -2.0, -2.0, -3.0, -3.0, -1.0, -2.0, -4.0, 7.0,-1.0, -1.0, -4.0, -3.0, -2.0, -2.0, -1.0, -2.0, -4.0],
                           ['S', 1.0, -1.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0, -2.0, -2.0, 0.0, -1.0, -2.0, -1.0, 4.0,1.0, -3.0, -2.0, -2.0, 0.0, 0.0, 0.0, -4.0],
                           ['T', 0.0, -1.0, 0.0, -1.0, -1.0, -1.0, -1.0, -2.0, -2.0, -1.0, -1.0, -1.0, -1.0, -2.0, -1.0,1.0, 5.0, -2.0, -2.0, 0.0, -1.0, -1.0, 0.0, -4.0],
                           ['W', -3.0, -3.0, -4.0, -4.0, -2.0, -2.0, -3.0, -2.0, -2.0, -3.0, -2.0, -3.0, -1.0, 1.0, -4.0,-3.0, -2.0, 11.0, 2.0, -3.0, -4.0, -3.0, -2.0, -4.0],
                           ['Y', -2.0, -2.0, -2.0, -3.0, -2.0, -1.0, -2.0, -3.0, 2.0, -1.0, -1.0, -2.0, -1.0, 3.0, -3.0,-2.0, -2.0, 2.0, 7.0, -1.0, -3.0, -2.0, -1.0, -4.0],
                           ['V', 0.0, -3.0, -3.0, -3.0, -1.0, -2.0, -2.0, -3.0, -3.0, 3.0, 1.0, -2.0, 1.0, -1.0, -2.0, -2.0,0.0, -3.0, -1.0, 4.0, -3.0, -2.0, -1.0, -4.0],
                           ['B', -2.0, -1.0, 3.0, 4.0, -3.0, 0.0, 1.0, -1.0, 0.0, -3.0, -4.0, 0.0, -3.0, -3.0, -2.0, 0.0,-1.0, -4.0, -3.0, -3.0, 4.0, 1.0, -1.0, -4.0],
                           ['Z', -1.0, 0.0, 0.0, 1.0, -3.0, 3.0, 4.0, -2.0, 0.0, -3.0, -3.0, 1.0, -1.0, -3.0, -1.0, 0.0,-1.0, -3.0, -2.0, -2.0, 1.0, 4.0, -1.0, -4.0],
                           ['X', 0.0, -1.0, -1.0, -1.0, -2.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -2.0,0.0, 0.0, -2.0, -1.0, -1.0, -1.0, -1.0, -1.0, -4.0],
                           ['-', -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0,-4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, 1.0],
                           ], dtype=object)
    return blosum_62_array


def create_dictionary_from_alignment(alignment_file):
    """Takes a fasta style alignment and makes a dictionary where the key is whatever signifier follows '>'
    and the value is the sequence with no spaces"""
    alignment = open(alignment_file, "r").read().split('>')
    sequence_dictionary = {sequence.split('\n')[0]: sequence.split('\n')[1].strip() for sequence in alignment if len(sequence) != 0}
    return sequence_dictionary


def shannons_entropy_for_residue_conservation(reference_sequence_fasta_name,alignment_file):
    """Takes a fasta style alignment creates a list of shannon's entropy values
    for each residue position with an amino acid in the reference sequence given based on blosum 62 matrix scores"""
    from numpy import where
    from math import log
    from collections import Counter
    # creating a dictionary of sequences where the keys are names signified in the fasta alignment alignment_file
    sequence_dictionary = create_dictionary_from_alignment(alignment_file)
    # Storing the sequence marked by reference_sequence_fasta_name
    reference_sequence = sequence_dictionary[reference_sequence_fasta_name]
    # making a list that gets rid of dashes in the sequence
    # and attaches the residue index in the entire alignment
    corrected_reference_sequence = [(residue_index,residue) for residue_index, residue in enumerate(reference_sequence) if residue != '-']
    # deleting the reference after storing for easier iteration later
    del sequence_dictionary[reference_sequence_fasta_name]
    # creating a list of the alignment names for iteration later
    sequence_list = list(sequence_dictionary.keys())
    shannons_entropy=[]
    # these sets of loops will work to calculate
    # and store entropy values for each residue position in the reference specified
    for residue_index,residue in corrected_reference_sequence:
        list_of_scores=[]
        # The residue letter code is stored in an array, and the corresponding row in the matrix is stored
        reference_matrix_index=where(blosum_62_matrix()[:,0]==reference_sequence[residue_index])[0]
        for name in sequence_list:
            # The comparison matrix column is stored, and used the row reference to
            # find the substitution score which is also stored,
            # this iterates for all comparison sequences
            comparison_matrix_index = where(blosum_62_matrix()[0,:] == sequence_dictionary[name][residue_index])[0]
            list_of_scores.append(int(blosum_62_matrix()[reference_matrix_index, comparison_matrix_index]))
        # The substitution scores for the reference amino acid are counted,
        # as well as all possible matrix values for the reference residue
        score_list_histogram=Counter(list_of_scores)
        number_of_possible_bins = len(Counter(list(blosum_62_matrix()[reference_matrix_index, 1:][0])))
        entropy = 0
        for score_count in score_list_histogram.values():
            # The entropy for each residue position is calculated and stored in a list
            probability = score_count/len(list_of_scores)
            entropy += -probability*log(probability, number_of_possible_bins)
        shannons_entropy.append(entropy)
    return shannons_entropy


def color_coding(protein_to_color,list_of_values, scale_tuple, color_scheme_tuple):
    """Takes a list_of_values and colors a residue of a pymol structure in the scale specified by
    scale_tuple:(min_value_in_scale,max_value_in_scale) and color_scheme_tuple:(low_color,high_color).
    This function has to be run in pymol to work and the protein that is being color coded has to already exist"""
    from pymol import cmd
    minimum_value=scale_tuple[0]
    maximum_value=scale_tuple[1]
    low_color=color_scheme_tuple[0]
    high_color=color_scheme_tuple[1]
    # using list comprehension to normalize the given values between 0 and 1 with the given scale
    normalized_values=[(original_value-minimum_value)/(maximum_value-minimum_value) for original_value in list_of_values]
    # assigning each normalized value to the b factor of each residue in the structure
    for residue_index, value in enumerate(normalized_values):
        cmd.alter('resi '+str(residue_index+1),'b='+str(value))
    # Lower values are assigned the first color listed, higher values latter color
    cmd.spectrum('b',low_color+'_'+high_color,protein_to_color)
    # spectrumbar green, blue, name = bar, head = (277.514, 232.761, 204.882), tail = (277.514, 252.761, 204.882), radius = 5


