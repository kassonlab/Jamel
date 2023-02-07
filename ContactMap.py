
def calc_residue_dist(residue_one, residue_two) :
    """Returns the C-alpha distance between two residues"""
    # Mostly copied from https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/protein_contact_map/
    from numpy import sqrt
    # Using distance formula to calculate residue distance
    diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord
    return sqrt(sum(diff_vector * diff_vector))
def residue_dist_matrix(pdb_file,chain_identifier,make_it_binary='Yes',distance_cutoff=7):
    """Takes a distance matrix of a protein and converts it to a binary matrix that uses 1 to signify a residue
    contact and 0 for no contact or Returns a matrix of C-alpha distances between residues in a protein chain or b"""
    # Mostly copied from https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/protein_contact_map/
    from numpy import array
    from Bio.PDB import PDBParser
    RESIDUE_NUMBER_CUTOFF=6
    # Calling the PDB structure with the parser
    STRUCTURE = PDBParser().get_structure(pdb_file, pdb_file)
    # Then selecting the chain
    PROTEIN = STRUCTURE[0][chain_identifier]
    if make_it_binary=='Yes':
        # Using list comprehension to create a nested list
        # then an array(protein length by protein length) of 1 to signifying a residue contact below distance_cutoff
        # and at least 6 residue separation or 0 for no contact
        return array([[1 if (calc_residue_dist(res_one, res_two)) <= distance_cutoff and (row-col) >= RESIDUE_NUMBER_CUTOFF else 0 for row, res_one in enumerate(PROTEIN)] for col, res_two in enumerate(PROTEIN)])
    else:
        # Using list comprehension to create a list then an array of each residue distance
        # to the rest of the chain and stacking those tuples to make a matrix
        return array([[calc_residue_dist(res_one,res_one) for res_one in PROTEIN] for res_one in PROTEIN])


def get_residue_contact_pairs(pdb_filename,chain_identifier):
    """Takes a pdb structure and return a nested list where each residue index in the list has a list with the python
    indexes for the residues it's in contact with.
    Residue indexes that have no contacts will have an empty list"""
    from numpy import where
    dist_matrix = residue_dist_matrix(pdb_filename,chain_identifier)
    x_axis,y_axis=list(where(dist_matrix==1)[0]),list(where(dist_matrix==1)[1])
    list_of_contact_pairs=[[] for rows in dist_matrix]
    for x, y in zip(x_axis, y_axis):
        list_of_contact_pairs[x].append(y)
    return list_of_contact_pairs


def correct_residue_position_for_alignment(pdb_file,chain_identifier,sequence_from_alignment):
    """Takes a nested list of residue positions
    and translates them into the residue position they are found in the alignment given"""
    contact_pairs = get_residue_contact_pairs(pdb_file,chain_identifier)
    sequence_indexing = [ind for ind, x in enumerate(sequence_from_alignment) if x != '-']
    index_dictionary = {real_index:alignment_index for real_index, alignment_index in enumerate(sequence_indexing)}
    updated_contact_map = [[index_dictionary[real_index] for real_index in contact_pairs[alignment_index]] for alignment_index, pairs in enumerate(contact_pairs)]
    return updated_contact_map


from ColorCoding import blosum_62_matrix,create_dictionary_from_alignment
from numpy import delete
from numpy import max as npmax
from numpy import min as npmin
# residuenumber,blosumscore,absolute value of confidence score difference,binary contact
blosum=blosum_62_matrix()
blosum=delete(blosum,0,axis=0)
blosum=delete(blosum,0,axis=1)
maxi, mini=npmax(blosum),npmin(blosum)
print(maxi, mini)
# normalize between zero an -1 and then flip the signs, also check distribution of blosum scores
# low + (high-low)*(x-minimum)/(maximum-minimmum)
normalize = lambda x: (-1+1*(x-mini)/(maxi-mini))
blosum=[list(normalize(y) for y in x) for x in blosum]
# 1 if both contact, 0 if no contact, -1 if only one contact

print(blosum)

def ContactOverlap(alignment_file, comparison, reference='6VSB_B'):
    from ColorCoding import create_dictionary_from_alignment
    from numpy import where
    sequence_dictionary=create_dictionary_from_alignment(alignment_file)
    reference_sequence,comparison_sequence=sequence_dictionary[reference],sequence_dictionary[comparison]
    #CP is Comparison Protein and RP is Reference Protein
    cp_updated_contact_map=correct_residue_position_for_alignment(f'3mer{comparison}.pdb','B',comparison_sequence)
    print(comparison_sequence)
    print(cp_updated_contact_map)
    rp_updated_contact_map = correct_residue_position_for_alignment(f'{reference}.pdb',"B", reference_sequence)
    rp_contact_map,cp_contact_map=[],[]
    j=0
    residue_comparison=[blosum[where(blosum_62_matrix()[:, 0] == res_ref)[0][0]-1][where(blosum_62_matrix()[0, :] == res_com)[0][0]-1]
                        for res_ref, res_com in zip(reference_sequence,comparison_sequence)]
    for residue in reference_sequence:
        if residue.isalpha():
            rp_contact_map.append(rp_updated_contact_map[j])
            j+=1
        else:
            rp_contact_map.append([])
    j=0
    for residue in comparison_sequence:
        if residue.isalpha():
            cp_contact_map.append(cp_updated_contact_map[j])
            j+=1
        else:
            cp_contact_map.append([])
    fraction_conserved=[]
    # print(list(zip(rp_contact_map,cp_contact_map)))

    for x, y in zip(rp_contact_map,cp_contact_map):
        if x==[] and y==[]: fraction_conserved.append(1); continue
        elif x == [] or y == []: fraction_conserved.append(0); continue
        x_set = set(); y_set = set()
        [x_set.update(num for num in range(contact-6,contact+7)) for contact in x]
        [y_set.update(num for num in range(contact - 6, contact + 7)) for contact in y]
        fraction_conserved.append((len([contact for contact in y if contact in x_set])+len([contact for contact in x if contact in y_set]))/(len(x)+len(y)))
    print(fraction_conserved)

ContactOverlap('SARS2wEverythingstable.aln','EidolonBat')
