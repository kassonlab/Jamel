from numpy import sqrt,array,where,delete,zeros,savetxt
from Bio import PDB,SeqIO
from ColorCoding import blosum_62_matrix,create_dictionary_from_alignment
from Analysis import rank_difference_table_to_dict
from numpy import max as npmax
from numpy import min as npmin
# TODO allow looking at contacts for specific position
# TODO combine contacts across all homomers
# TODO deprecate other contact finding method with kdtree
from Bio import PDB
from scipy.spatial import KDTree

def get_coordinates_dict_per_chain(pdb_file):
    STRUCTURE = PDB.PDBParser().get_structure(pdb_file, pdb_file)
    # Then selecting the chain
    PROTEIN = STRUCTURE[0]
    list_of_chain_coordinates = {chain.id: [tuple(res["CA"].coord) for res in chain] for chain in PROTEIN}
    return list_of_chain_coordinates


def create_spatial_index(list_coords):
    kdtree = KDTree(list_coords)
    return kdtree


def find_neighbors(kdtree, list_of_coords, residue_index, distance_cutoff=7.0, residue_dist_cutoff=6):
    neighbors = kdtree.query_ball_point(list_of_coords[residue_index], distance_cutoff)
    return {res for res in neighbors if res - residue_index >= residue_dist_cutoff}


def assign_index_to_chain(index, contact_range):
    for chain, ranges in contact_range.items():
        if index in ranges:
            corrected_index=contact_range[chain].index(index)
            return chain,corrected_index


def intra_protein_contacts(pdb_file, chain_id):
    coord_dict = get_coordinates_dict_per_chain(pdb_file)
    contact_range = {}
    range_placeholder = 0
    total_coords = []
    for chain, coords in coord_dict.items():
        contact_range[chain] = range(range_placeholder, range_placeholder + len(coords))
        total_coords += coords
        range_placeholder += len(coords)
    kd_tree = create_spatial_index(total_coords)
    intra_contacts = []
    for index,kd_index in enumerate(contact_range[chain_id]):
        neighbors = find_neighbors(kd_tree, total_coords, kd_index, residue_dist_cutoff=0)
        chain_assigned_neighbors = []
        for neighbor in neighbors:
            chain_for_index=assign_index_to_chain(neighbor, contact_range)
            if chain_for_index[0]!=chain_id:
                chain_assigned_neighbors.append(f'Ch{chain_for_index[0]}:{chain_for_index[1]}')
        intra_contacts.append(chain_assigned_neighbors)
    return intra_contacts

def calc_residue_dist(residue_one, residue_two) :
    """Returns the C-alpha distance between two residues"""
    # Mostly copied from https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/protein_contact_map/

    # Using distance formula to calculate residue distance
    diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord
    return sqrt(sum(diff_vector * diff_vector))
def residue_dist_matrix(pdb_file,chain_identifier,make_it_binary='Yes',distance_cutoff=7):
    """Takes a distance matrix of a protein and converts it to a binary matrix that uses 1 to signify a residue
    contact and 0 for no contact or Returns a matrix of C-alpha distances between residues in a protein chain or b"""
    # Mostly copied from https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/protein_contact_map/

    RESIDUE_NUMBER_CUTOFF=6
    # Calling the PDB structure with the parser
    STRUCTURE = PDB.PDBParser().get_structure(pdb_file, pdb_file)
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
        return array([[calc_residue_dist(res_one,res_two) for res_one in PROTEIN] for res_two in PROTEIN])


def get_residue_contact_pairs(pdb_filename,chain_identifier):
    """Takes a pdb structure and return a nested list where each residue index in the list has a list with the python
    indexes for the residues it's in contact with.
    Residue indexes that have no contacts will have an empty list"""

    dist_matrix = residue_dist_matrix(pdb_filename,chain_identifier)
    x_axis,y_axis=list(where(dist_matrix==1)[0]),list(where(dist_matrix==1)[1])
    list_of_contact_pairs=[[] for rows in dist_matrix]
    for x, y in zip(x_axis, y_axis):
        list_of_contact_pairs[x].append(y)
    return list_of_contact_pairs


def correct_alignment_for_residue_position(alignment_file,alignment_label,alignment_position):
    sequence_dict=create_dictionary_from_alignment(alignment_file)
    sequence_from_alignment=sequence_dict[alignment_label]
    sequence_indexing = [ind for ind, x in enumerate(sequence_from_alignment) if x != '-']
    return sequence_indexing.index(alignment_position)


def get_residue_at_native_position(alignment_file,alignment_label,alignment_index=None,residue_index=None):
    sequence_dict=create_dictionary_from_alignment(alignment_file)
    sequence_from_alignment=sequence_dict[alignment_label]
    if alignment_index:
        return sequence_from_alignment[alignment_index]
    elif residue_index:
        sequence = [x for ind, x in enumerate(sequence_from_alignment) if x != '-']
        return sequence[residue_index]


def correct_residue_position_for_alignment(pdb_file,chain_identifier,alignment_file,alignment_label):
    """This can only be used for the pdbs that have the same sequence as the sequence found in the alignment. Takes a nested list of residue indexes
    and translates them into the alignment indexes they are found in the alignment given"""
    sequence_from_alignment=create_dictionary_from_alignment(alignment_file)[alignment_label]
    contact_pairs = get_residue_contact_pairs(pdb_file,chain_identifier)
    sequence_indexing = [ind for ind, x in enumerate(sequence_from_alignment) if x != '-']
    index_dictionary = {real_index:alignment_index for real_index, alignment_index in enumerate(sequence_indexing)}
    updated_contact_map = [[index_dictionary[real_index] for real_index in contact_pairs[alignment_index]] for alignment_index, pairs in enumerate(contact_pairs)]
    return updated_contact_map
def correct_native_for_chimera_index(alignment_file,label,native_residue_index,chimera_seq):
    sequence_from_alignment = create_dictionary_from_alignment(alignment_file)[label]
    chunk_to_find=sequence_from_alignment.replace('-','')[native_residue_index:native_residue_index+4]
    return chimera_seq.find(chunk_to_find)

def get_sequence_from_pdb(pdb_file, chain_id):
    # Iterate over all records in the PDB file
    for record in SeqIO.parse(pdb_file, "pdb-atom"):
        # Check if the record corresponds to the desired chain
        if record.annotations["chain"] == chain_id:
            # Return the sequence
            return str(record.seq)
    # Extract the sequence from the chain


def compare_contacts(alignment_file,native_pdb,chimera_pdb,chain_identifier,label,alignment_index,rank_difference_file):
    try:
        residue=get_residue_at_native_position(alignment_file,protein_label,alignment_index)
        native_index = correct_alignment_for_residue_position(alignment_file, label, alignment_index)
    except:
        print('doesnt exist')
        return [label, alignment_index + 1, [], [], [], []]
    native_contacts = get_residue_contact_pairs(native_pdb, chain_identifier)[native_index]
    native_seq = get_sequence_from_pdb(native_pdb, chain_identifier)
    if native_contacts:
        native_contacts_ids = [f'{native_seq[x]}{x + 1}' for x in native_contacts]
    else:
        native_contacts_ids = []
    chimera_seq = get_sequence_from_pdb(chimera_pdb, chain_identifier)
    chi_index = correct_native_for_chimera_index(alignment_file, protein_label, native_index, chimera_seq)
    residue_ids = f'{residue}{native_index + 1},{residue}{chi_index + 1}'
    chimera_contacts = get_residue_contact_pairs(chimera_pdb, chain_identifier)[chi_index]
    if chimera_contacts:
        chimera_contacts_ids = [f'{chimera_seq[x]}{x + 1}' for x in chimera_contacts]
    else:
        chimera_contacts_ids = []
    rank_difference = rank_difference_table_to_dict(rank_difference_file)[label]
    return [label, alignment_index + 1, residue_ids, rank_difference, native_contacts_ids, chimera_contacts_ids]


with open("/gpfs/gpfs0/scratch/jws6pq/Notebook/Overall/List_of_coronaviruses", 'r') as loc:
    loc = loc.readlines()
indexes=['Human229E',
'WigeonHKU20',
'Wencheng',
'SorexT14',
'EidolonBat',
'BATGCCDC1',
'BatBGR',
'Shandong',
'BetaCoronaSC2018',
'Mystacina','RabbitHKU14','BatAlpha','MoorHKU21','FalconHKU27','WIV16']
alignment_positions=[1262,1263,1264,1265,1266,1267,1268]
aln='/gpfs/gpfs0/scratch/jws6pq/Notebook/PDB/CoronavirusMSA.aln'
comparison_matrix=zeros((len(alignment_positions)*len(indexes),6),dtype=object)
y=0
for postion in alignment_positions:
    rank_change = f'/gpfs/gpfs0/scratch/jws6pq/Notebook/PDB/Rank_change_{postion}.tsv'
    for label in indexes:
        protein_label=label
        native=f'/gpfs/gpfs0/scratch/jws6pq/Notebook/PDB/3mer{protein_label}.pdb'
        chi=f'/gpfs/gpfs0/scratch/jws6pq/Notebook/PDB/3merSARS2w{protein_label}S1.pdb'
        comparison=compare_contacts(aln, native, chi, 'B', protein_label, postion - 1, rank_change)
        print(comparison)
        comparison_matrix[y]=comparison
        y+=1
savetxt(f'/gpfs/gpfs0/scratch/jws6pq/Notebook/PDB/rank_comparison.csv', comparison_matrix, fmt='%s-%s-%s-%s-%s-%s')
# residuenumber,blosumscore,absolute value of confidence score difference,binary contact
# blosum=blosum_62_matrix()
# blosum=delete(blosum,0,axis=0)
# blosum=delete(blosum,0,axis=1)
# maxi, mini=npmax(blosum),npmin(blosum)
# normalize between zero an -1 and then flip the signs, also check distribution of blosum scores
# low + (high-low)*(x-minimum)/(maximum-minimmum)
# normalize = lambda x: (-1+1*(x-mini)/(maxi-mini))
# blosum=[list(normalize(y) for y in x) for x in blosum]
# 1 if both contact, 0 if no contact, -1 if only one contact


def ContactOverlap(alignment_file, comparison, reference='6VSB_B'):
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

