
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
        # Using list comprehension to create a list then an array of each residue distance to the rest of the chain and stacking those tuples to make a matrix
        return array([[calc_residue_dist(x,y) for x in PROTEIN] for y in PROTEIN])


def get_residue_contact_pairs(pdb_filename,chain_identifier):
    """Takes a pdb structure and return a nested list where each residue index in the list has a list with the python
    indexes for the residues it's in contact with.
    Residue indexes that have no contacts will have an empty list"""
    from numpy import where
    dist_matrix = residue_dist_matrix(pdb_filename,chain_identifier)
    x_axis,y_axis=list(where(dist_matrix==1)[0]),list(where(dist_matrix==1)[1])
    list_of_contact_pairs=[[] for x in dist_matrix]
    for x, y in zip(x_axis, y_axis):
        list_of_contact_pairs[x].append(y)
    return list_of_contact_pairs
print(get_residue_contact_pairs('6VSB_B.pdb','B'))


def correct_residue_position_for_alignment(pdb_file,chain_identifier,alignment):
    """Takes an alignment and creates a nested list"""
    contact_map = get_residue_contact_pairs(pdb_file,chain_identifier)
    sequence_indexing = [ind for ind, x in enumerate(alignment) if x != '-']
    residue_position_dictionary = {indx: indy for indx, indy in enumerate(sequence_indexing)}
    updated_contact_map = [[residue_position_dictionary[y] for y in contact_map[ind]] for ind, x in enumerate(contact_map)]
    return updated_contact_map
correct_residue_position_for_alignment()

def ContactOverlap(Alignmentfile,comparison,reference='6vsb_B'):
    # Alignment in FASTA format. Make sure your benchmark sequence is first
    Sequences = open(Alignmentfile, "r").read().split('>')
    SequenceDictionary={sequence.split('\n')[0]:sequence.split('\n')[1].strip() for sequence in Sequences if len(sequence)!=0}
    ReferenceSequence,ComparisonSequence=SequenceDictionary[reference],SequenceDictionary[comparison]
    #CP is Comparison Protein and RP is Reference Protein
    CPUpdatedContactMap=CorrectResiduePositionforAlignment(comparison,ComparisonSequence)
    RPUpdatedContactMap = CorrectResiduePositionforAlignment(reference, ReferenceSequence)
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
