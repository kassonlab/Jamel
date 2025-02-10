from Bio import SeqIO
from numpy import array
from math import log
from collections import Counter
from Chimeragenesis.AccessiontoAlignment import create_dictionary_from_alignment
from Chimeragenesis.ChimeraGenerator import fasta_creation
def blosum_62_matrix():
    """Simply returns the blosum 62 matrix"""

    blosum_62_array = array([['', 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T',
                              'W', 'Y', 'V', 'B', 'Z', 'X', '-'],
                             ['A', 4.0, -1.0, -2.0, -2.0, 0.0, -1.0, -1.0, 0.0, -2.0, -1.0, -1.0, -1.0, -1.0, -2.0,
                              -1.0, 1.0, 0.0, -3.0, -2.0, 0.0, -2.0, -1.0, 0.0, -4.0],
                             ['R', -1.0, 5.0, 0.0, -2.0, -3.0, 1.0, 0.0, -2.0, 0.0, -3.0, -2.0, 2.0, -1.0, -3.0, -2.0,
                              -1.0, -1.0, -3.0, -2.0, -3.0, -1.0, 0.0, -1.0, -4.0],
                             ['N', -2.0, 0.0, 6.0, 1.0, -3.0, 0.0, 0.0, 0.0, 1.0, -3.0, -3.0, 0.0, -2.0, -3.0, -2.0,
                              1.0, 0.0, -4.0, -2.0, -3.0, 3.0, 0.0, -1.0, -4.0],
                             ['D', -2.0, -2.0, 1.0, 6.0, -3.0, 0.0, 2.0, -1.0, -1.0, -3.0, -4.0, -1.0, -3.0, -3.0, -1.0,
                              0.0, -1.0, -4.0, -3.0, -3.0, 4.0, 1.0, -1.0, -4.0],
                             ['C', 0.0, -3.0, -3.0, -3.0, 9.0, -3.0, -4.0, -3.0, -3.0, -1.0, -1.0, -3.0, -1.0, -2.0,
                              -3.0, -1.0, -1.0, -2.0, -2.0, -1.0, -3.0, -3.0, -2.0, -4.0],
                             ['Q', -1.0, 1.0, 0.0, 0.0, -3.0, 5.0, 2.0, -2.0, 0.0, -3.0, -2.0, 1.0, 0.0, -3.0, -1.0,
                              0.0, -1.0, -2.0, -1.0, -2.0, 0.0, 3.0, -1.0, -4.0],
                             ['E', -1.0, 0.0, 0.0, 2.0, -4.0, 2.0, 5.0, -2.0, 0.0, -3.0, -3.0, 1.0, -2.0, -3.0, -1.0,
                              0.0, -1.0, -3.0, -2.0, -2.0, 1.0, 4.0, -1.0, -4.0],
                             ['G', 0.0, -2.0, 0.0, -1.0, -3.0, -2.0, -2.0, 6.0, -2.0, -4.0, -4.0, -2.0, -3.0, -3.0,
                              -2.0, 0.0, -2.0, -2.0, -3.0, -3.0, -1.0, -2.0, -1.0, -4.0],
                             ['H', -2.0, 0.0, 1.0, -1.0, -3.0, 0.0, 0.0, -2.0, 8.0, -3.0, -3.0, -1.0, -2.0, -1.0, -2.0,
                              -1.0, -2.0, -2.0, 2.0, -3.0, 0.0, 0.0, -1.0, -4.0],
                             ['I', -1.0, -3.0, -3.0, -3.0, -1.0, -3.0, -3.0, -4.0, -3.0, 4.0, 2.0, -3.0, 1.0, 0.0, -3.0,
                              -2.0, -1.0, -3.0, -1.0, 3.0, -3.0, -3.0, -1.0, -4.0],
                             ['L', -1.0, -2.0, -3.0, -4.0, -1.0, -2.0, -3.0, -4.0, -3.0, 2.0, 4.0, -2.0, 2.0, 0.0, -3.0,
                              -2.0, -1.0, -2.0, -1.0, 1.0, -4.0, -3.0, -1.0, -4.0],
                             ['K', -1.0, 2.0, 0.0, -1.0, -3.0, 1.0, 1.0, -2.0, -1.0, -3.0, -2.0, 5.0, -1.0, -3.0, -1.0,
                              0.0, -1.0, -3.0, -2.0, -2.0, 0.0, 1.0, -1.0, -4.0],
                             ['M', -1.0, -1.0, -2.0, -3.0, -1.0, 0.0, -2.0, -3.0, -2.0, 1.0, 2.0, -1.0, 5.0, 0.0, -2.0,
                              -1.0, -1.0, -1.0, -1.0, 1.0, -3.0, -1.0, -1.0, -4.0],
                             ['F', -2.0, -3.0, -3.0, -3.0, -2.0, -3.0, -3.0, -3.0, -1.0, 0.0, 0.0, -3.0, 0.0, 6.0, -4.0,
                              -2.0, -2.0, 1.0, 3.0, -1.0, -3.0, -3.0, -1.0, -4.0],
                             ['P', -1.0, -2.0, -2.0, -1.0, -3.0, -1.0, -1.0, -2.0, -2.0, -3.0, -3.0, -1.0, -2.0, -4.0,
                              7.0, -1.0, -1.0, -4.0, -3.0, -2.0, -2.0, -1.0, -2.0, -4.0],
                             ['S', 1.0, -1.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0, -2.0, -2.0, 0.0, -1.0, -2.0, -1.0,
                              4.0, 1.0, -3.0, -2.0, -2.0, 0.0, 0.0, 0.0, -4.0],
                             ['T', 0.0, -1.0, 0.0, -1.0, -1.0, -1.0, -1.0, -2.0, -2.0, -1.0, -1.0, -1.0, -1.0, -2.0,
                              -1.0, 1.0, 5.0, -2.0, -2.0, 0.0, -1.0, -1.0, 0.0, -4.0],
                             ['W', -3.0, -3.0, -4.0, -4.0, -2.0, -2.0, -3.0, -2.0, -2.0, -3.0, -2.0, -3.0, -1.0, 1.0,
                              -4.0, -3.0, -2.0, 11.0, 2.0, -3.0, -4.0, -3.0, -2.0, -4.0],
                             ['Y', -2.0, -2.0, -2.0, -3.0, -2.0, -1.0, -2.0, -3.0, 2.0, -1.0, -1.0, -2.0, -1.0, 3.0,
                              -3.0, -2.0, -2.0, 2.0, 7.0, -1.0, -3.0, -2.0, -1.0, -4.0],
                             ['V', 0.0, -3.0, -3.0, -3.0, -1.0, -2.0, -2.0, -3.0, -3.0, 3.0, 1.0, -2.0, 1.0, -1.0, -2.0,
                              -2.0, 0.0, -3.0, -1.0, 4.0, -3.0, -2.0, -1.0, -4.0],
                             ['B', -2.0, -1.0, 3.0, 4.0, -3.0, 0.0, 1.0, -1.0, 0.0, -3.0, -4.0, 0.0, -3.0, -3.0, -2.0,
                              0.0, -1.0, -4.0, -3.0, -3.0, 4.0, 1.0, -1.0, -4.0],
                             ['Z', -1.0, 0.0, 0.0, 1.0, -3.0, 3.0, 4.0, -2.0, 0.0, -3.0, -3.0, 1.0, -1.0, -3.0, -1.0,
                              0.0, -1.0, -3.0, -2.0, -2.0, 1.0, 4.0, -1.0, -4.0],
                             ['X', 0.0, -1.0, -1.0, -1.0, -2.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,
                              -2.0, 0.0, 0.0, -2.0, -1.0, -1.0, -1.0, -1.0, -1.0, -4.0],
                             ['-', -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0,
                              -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, 1.0],
                             ], dtype=object)
    return blosum_62_array





def shannons_entropy_for_residue_conservation(reference_label, alignment_file):
    """Takes a fasta style alignment creates a list of shannon's entropy values
    for each residue position with an amino acid in the reference sequence given based on blosum 62 matrix scores"""

    # creating a dictionary of sequences where the keys are names signified in the fasta alignment alignment_file
    sequence_dictionary = create_dictionary_from_alignment(alignment_file)
    # Storing the sequence marked by reference_sequence_fasta_name
    reference_sequence = sequence_dictionary[reference_label]
    # making a list that gets rid of dashes in the sequence
    # and attaches the residue index in the entire alignment
    alignment_indexed_reference = [(alignment_index, residue) for alignment_index, residue in enumerate(reference_sequence)
                                    if residue != '-']
    # deleting the reference after storing for easier iteration later
    del sequence_dictionary[reference_label]
    # creating a list of the alignment names for iteration later
    shannons_entropy = []
    # these sets of loops will work to calculate
    # and store entropy values for each residue position in the reference specified
    blosum_column=list(blosum_62_matrix()[:,0])
    blosum_row=list(blosum_62_matrix()[0,:])
    for alignment_index, residue in alignment_indexed_reference:
        list_of_scores = []
        # The residue letter code is stored in an array, and the corresponding row in the matrix is stored
        reference_blosum_row = blosum_row.index(reference_sequence[alignment_index])
        for sequence in sequence_dictionary.values():
            # The comparison matrix column is stored, and used the row reference to
            # find the substitution score which is also stored,
            # this iterates for all comparison sequences
            comparison_blosum_column = blosum_column.index(sequence[alignment_index])
            list_of_scores.append(int(blosum_62_matrix()[reference_blosum_row, comparison_blosum_column]))
        # The substitution scores for the reference amino acid are counted,
        # as well as all possible matrix values for the reference residue
        score_list_histogram = Counter(list_of_scores)
        number_of_possible_bins = len(Counter(list(blosum_62_matrix()[reference_blosum_row, 1:])))
        entropy = 0
        for score_count in score_list_histogram.values():
            # The entropy for each residue position is calculated and stored in a list
            probability = score_count / len(list_of_scores)
            entropy += -probability * log(probability, number_of_possible_bins)
        shannons_entropy.append((alignment_index,reference_sequence[alignment_index],entropy))
    return shannons_entropy
def color_coding(pdb_name, list_of_values, scale_tuple, color_scheme_tuple):
    """Takes a list_of_values and colors a residue of a pymol structure in the scale specified by
    scale_tuple:(min_value_in_scale,max_value_in_scale) and color_scheme_tuple:(low_color,high_color).
    This function has to be run in pymol to work and the protein that is being color coded has to already exist"""
    from pymol import cmd
    minimum_value = scale_tuple[0]
    maximum_value = scale_tuple[1]
    low_color = color_scheme_tuple[0]
    high_color = color_scheme_tuple[1]
    # using list comprehension to normalize the given values between 0 and 1 with the given scale
    normalized_values = [(original_value - minimum_value) / (maximum_value - minimum_value) for original_value in
                         list_of_values]
    # assigning each normalized value to the b factor of each residue in the structure
    for residue_index, value in enumerate(normalized_values):
        cmd.alter('resi ' + str(residue_index + 1), 'b=' + str(value))
    # Lower values are assigned the first color listed, higher values latter color
    cmd.spectrum('b', low_color + '_' + high_color, pdb_name)
    # spectrumbar green, blue, name = bar, head = (277.514, 232.761, 204.882), tail = (277.514, 252.761, 204.882), radius = 5


def plddt_color_coding(plddt_file, pdb_name, color_scheme_tuple):
    """Takes a list_of_values and colors a residue of a pymol structure in the scale specified by
    and color_scheme_tuple:(low_color,high_color) according to the corresponding plddt value given in plddt_file.
    This function has to be run in pymol to work and the protein that is being color coded has to already exist"""
    list_of_values = [float(score) for score in open(plddt_file, 'r').readlines()]
    color_coding(pdb_name, list_of_values, (0, 100), color_scheme_tuple)
def pdb_to_fasta(pdb,new_fasta):
    with open(pdb, 'r') as pdb_file:
        fasta_list=[]
        for record in SeqIO.parse(pdb_file, 'pdb-atom'):
            fasta_list.append((record.seq,1,record.id))
    fasta_creation(new_fasta,fasta_list)
