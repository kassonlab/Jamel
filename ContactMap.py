
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


from ColorCoding import blosum_62_matrix
from numpy import delete
from numpy import max as npmax
from numpy import min as npmin
# residuenumber,blosumscore,absolute value of confidence score difference,binary contact
blosum=blosum_62_matrix()
blosum=delete(blosum,0,axis=0)
blosum=delete(blosum,0,axis=1)
maxi, mini=npmax((blosum)),npmin(blosum)
print(maxi,mini)
lambda x:(2*((x-mini)/(maxi-mini))-1)*-1
print(blosum)

def ContactOverlap(Alignmentfile,comparison,reference='6vsb_B'):
    Sequences = open(Alignmentfile, "r").read().split('>')
    SequenceDictionary={sequence.split('\n')[0]:sequence.split('\n')[1].strip() for sequence in Sequences if len(sequence)!=0}
    ReferenceSequence,ComparisonSequence=SequenceDictionary[reference],SequenceDictionary[comparison]
    #CP is Comparison Protein and RP is Reference Protein
    CPUpdatedContactMap=correct_residue_position_for_alignment(comparison,ComparisonSequence)
    RPUpdatedContactMap = correct_residue_position_for_alignment(reference, ReferenceSequence)
    ReferenceContactMap,ComparisonContactMap=[],[]
    j=0
    #Should this be a function?
    RPContactCount,CPContactCount=0,0
    for x in ReferenceSequence:
        if x.isalpha():
            ReferenceContactMap.append([x]+RPUpdatedContactMap[j])
            RPContactCount+=len(RPUpdatedContactMap[j])+1
            j+=1
        else:
            ReferenceContactMap.append([x])
    j=0
    for x in ComparisonSequence:
        if x.isalpha():
            ComparisonContactMap.append([x]+CPUpdatedContactMap[j])
            CPContactCount+=len(CPUpdatedContactMap[j])+1
            j+=1
        else:
            ComparisonContactMap.append([x])
    TotalContacts = RPContactCount+CPContactCount
    for x,y in zip(ReferenceContactMap,ComparisonContactMap):
        for w in x:
            if w not in y and w!='-' and len(x)>1 or w not in y and w!='-' and len(y)>1:
                TotalContacts+=-1
        for v in y:
            if v not in x and v!='-' and len(x)>1 or v not in x and v!='-' and len(y)>1:
                TotalContacts += -1
    # system('/scratch/jws6pq/EMBOSS-6.6.0/emboss/needle -sprotein -gapopen 10 -gapextend 0.5 -outfile /gpfs/gpfs0/scratch/jws6pq/Notebook/Emboss/Full' + key + '.emboss -asequence /gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/' + key + '.fasta -bsequence /gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/SARS2.fasta')
    EmbossScore = open('/gpfs/gpfs0/scratch/jws6pq/Notebook/Emboss/Full' + comparison + '.emboss', 'r').readlines()[25].split()[-1]
    return comparison,TotalContacts,EmbossScore.replace('(','').replace(')','').replace('%','')
#Do i consider all the times where there are residues beyond SARS???????
#
#
