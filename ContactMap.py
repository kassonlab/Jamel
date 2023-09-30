from numpy import sqrt, array, where,  zeros
from Bio import PDB, SeqIO
from ColorCoding import blosum_62_matrix
from Analysis import convert_array_to_file
# TODO allow looking at contacts for specific position for intra
# TODO combine contacts across all homomers
from scipy.spatial import KDTree
from AccessiontoAlignment import get_alignment_indexing,get_alignment_indexing_w_dashes, alignment_finder,create_dictionary_from_alignment,no_gap_sequence_from_alignment


def get_coordinates_dict_per_chain(pdb_file):
    STRUCTURE = PDB.PDBParser().get_structure(pdb_file, pdb_file)
    # Then selecting the chain
    PROTEIN = STRUCTURE[0]
    list_of_chain_coordinates = {
        chain.id: [tuple(res["CB"].coord) if res.get_resname() != 'GLY' else tuple(res["CA"].coord) for res in chain if
                   len(chain) > 10] for chain in PROTEIN}

    return list_of_chain_coordinates


def create_spatial_index(list_coords):
    kdtree = KDTree(list_coords)
    return kdtree


def find_neighbors(kdtree, list_of_coords, residue_index, distance_cutoff=8.0, residue_dist_cutoff=6):
    neighbors = kdtree.query_ball_point(list_of_coords[residue_index], distance_cutoff)
    return sorted(tuple(res for res in neighbors if abs(res - residue_index) >= residue_dist_cutoff))


def assign_index_to_chain(index, contact_range):
    for chain, ranges in contact_range.items():
        if index in ranges:
            corrected_index = contact_range[chain].index(index)
            return chain, corrected_index


# TODO be able to add residue id
def get_inter_protein_contacts(pdb_file, chain_id, index_to_position=False):
    coord_dict = get_coordinates_dict_per_chain(pdb_file)
    chain_indices_dict = {}
    range_placeholder = 0
    total_coords = []
    for chain, coords in coord_dict.items():
        chain_indices_dict[chain] = range(range_placeholder, range_placeholder + len(coords))
        total_coords += coords
        range_placeholder += len(coords)
    kd_tree = create_spatial_index(total_coords)
    inter_contacts = []
    for index, kd_index in enumerate(chain_indices_dict[chain_id]):
        neighbors = find_neighbors(kd_tree, total_coords, kd_index, residue_dist_cutoff=0)
        chain_assigned_neighbors = []
        for neighbor in neighbors:
            chain_for_index = assign_index_to_chain(neighbor, chain_indices_dict)
            if chain_for_index[0] != chain_id:
                chain_assigned_neighbors.append(f'{chain_for_index[0]}:{chain_for_index[1] + bool(index_to_position)}')
        inter_contacts.append(chain_assigned_neighbors)
    return inter_contacts


def calc_residue_dist(residue_one, residue_two):
    """Returns the C-alpha distance between two residues"""
    # Mostly copied from https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/protein_contact_map/

    # Using distance formula to calculate residue distance
    diff_vector = residue_one["CB"].coord - residue_two["CB"].coord
    return sqrt(sum(diff_vector * diff_vector))


def residue_dist_matrix(pdb_file, chain_id, make_it_binary='Yes', distance_cutoff=8):
    """Takes a distance matrix of a protein and converts it to a binary matrix that uses 1 to signify a residue
    contact and 0 for no contact or Returns a matrix of C-alpha distances between residues in a protein chain or b"""

    if make_it_binary == 'Yes':
        # Using list comprehension to create a nested list
        # then an array(protein length by protein length) of 1 to signifying a residue contact below distance_cutoff
        # and at least 6 residue separation or 0 for no contact
        protein_coords = get_coordinates_dict_per_chain(pdb_file)[chain_id]
        contact_map = zeros((len(protein_coords), len(protein_coords)))
        contact_pairs = get_intra_residue_contact_pairs(pdb_file, chain_id)
        for index, pairs in enumerate(contact_pairs):
            for pair in pairs:
                contact_map[pair, index] = 1
                contact_map[index, pair] = 1
        return contact_map
    else:
        STRUCTURE = PDB.PDBParser().get_structure(pdb_file, pdb_file)
        # Then selecting the chain
        PROTEIN = STRUCTURE[0][chain_id]
        # Using list comprehension to create a list then an array of each residue distance
        # to the rest of the chain and stacking those tuples to make a matrix
        return array([[calc_residue_dist(res_one, res_two) for res_one in PROTEIN] for res_two in PROTEIN])


def get_intra_residue_contact_pairs(pdb_file, chain_id):
    """Takes a pdb structure and return a nested list where each residue index in the list has a list with the python
    indexes for the residues it's in contact with.
    Residue indexes that have no contacts will have an empty list"""

    protein_coords = get_coordinates_dict_per_chain(pdb_file)[chain_id]
    kd_tree = create_spatial_index(protein_coords)
    list_of_contact_pairs = [find_neighbors(kd_tree, protein_coords, index) for index, coord in
                             enumerate(protein_coords)]
    return list_of_contact_pairs


def correct_alignment_for_residue_position(alignment_file, alignment_label, alignment_index):
    sequence_dict = create_dictionary_from_alignment(alignment_file)
    sequence_from_alignment = sequence_dict[alignment_label]
    sequence_indexing = get_alignment_indexing(sequence_from_alignment)
    return sequence_indexing.index(alignment_index)


def correct_alignment_for_chimera_residue_position(aln_file, ref_label, aln_label, aln_index, spliced_in_seq):
    sequence_from_aln = create_dictionary_from_alignment(aln_file)[aln_label]
    ref_aln = create_dictionary_from_alignment(aln_file)[ref_label]
    boundaries = alignment_finder(aln_file, spliced_in_seq, aln_label, ref_label)
    com_boundaries = boundaries[1]
    com_boundaries = boundaries[2]
    comparison_aln_indexing = get_alignment_indexing(sequence_from_aln)[boundaries[1]]
    ref_aln_indexing = get_alignment_indexing(ref_aln)
    # nonspliced=
    ref_aln = create_dictionary_from_alignment(aln_file)[ref_label]


def correct_residue_index_for_alignment(alignment_file, residue_index, protein_label):
    aligned_seq = create_dictionary_from_alignment(alignment_file)[protein_label]
    print(get_alignment_indexing(aligned_seq)[residue_index])


def get_residue_at_native_position(alignment_file, alignment_label, alignment_index=None, residue_index=None):
    sequence_dict = create_dictionary_from_alignment(alignment_file)
    sequence_from_alignment = sequence_dict[alignment_label]
    if alignment_index:
        return sequence_from_alignment[alignment_index]
    elif residue_index:
        sequence = get_alignment_indexing(sequence_from_alignment)
        return sequence[residue_index]


def get_individual_intra_contacts(pdb_file, chain_id, index):
    coords = get_coordinates_dict_per_chain(pdb_file)[chain_id]
    kd_tree = create_spatial_index(coords)
    contacts = find_neighbors(kd_tree, coords, index)
    return contacts


def get_individual_inter_contacts(pdb_file, chain_id, index, index_to_position=False):
    coord_dict = get_coordinates_dict_per_chain(pdb_file)
    contact_range = {}
    range_placeholder = 0
    total_coords = []
    for chain, coords in coord_dict.items():
        contact_range[chain] = range(range_placeholder, range_placeholder + len(coords))
        total_coords += coords
        range_placeholder += len(coords)
    kd_tree = create_spatial_index(total_coords)
    inter_contacts = []
    neighbors = kd_tree.query_ball_point(coord_dict[chain_id][index], 8)
    for neighbor in neighbors:
        chain_for_index = assign_index_to_chain(neighbor, contact_range)
        if chain_for_index[0] != chain_id:
            inter_contacts.append(f'{chain_for_index[0]}:{chain_for_index[1] + bool(index_to_position)}')
    return inter_contacts


def correct_contact_position_for_alignment(pdb_file, chain_identifier, alignment_file, alignment_label):
    """This can only be used for the pdbs that have the same sequence as the sequence found in the alignment. Takes a nested list of residue indexes
    and translates them into the alignment indexes they are found in the alignment given"""
    sequence_from_alignment = create_dictionary_from_alignment(alignment_file)[alignment_label]
    contact_pairs = get_intra_residue_contact_pairs(pdb_file, chain_identifier)
    sequence_indexing = [ind for ind, x in enumerate(sequence_from_alignment) if x != '-']
    index_dictionary = {real_index: alignment_index for real_index, alignment_index in enumerate(sequence_indexing)}
    updated_contact_map = [[index_dictionary[real_index] for real_index in contact_pairs[alignment_index]] for
                           alignment_index, pairs in enumerate(contact_pairs)]
    return updated_contact_map


def correct_alignment_for_chimera_index(alignment_file, ref_label, aln_index, comparison_label,sequence_of_interest):
    sequence_dictionary = create_dictionary_from_alignment(alignment_file)
    reference_aln = sequence_dictionary[ref_label]
    comparison_aln = sequence_dictionary[comparison_label]
    reference_alignment_indexing = get_alignment_indexing(reference_aln)
    no_gap_reference_sequence = no_gap_sequence_from_alignment(reference_aln)
    alignment_reference_start = reference_alignment_indexing[no_gap_reference_sequence.find(sequence_of_interest)]
    alignment_reference_end = reference_alignment_indexing[no_gap_reference_sequence.find(sequence_of_interest) + len(
        sequence_of_interest) - 1] + 1
    indexed_ref_aln=get_alignment_indexing_w_dashes(reference_aln)
    indexed_com_aln=get_alignment_indexing_w_dashes(comparison_aln)
    indexed_ref_aln[alignment_reference_start:alignment_reference_end]=indexed_com_aln[alignment_reference_start:alignment_reference_end]
    ref_aln_list=list(reference_aln)
    com_aln_list=list(comparison_aln)
    ref_aln_list[alignment_reference_start:alignment_reference_end] = com_aln_list[
                                                                         alignment_reference_start:alignment_reference_end]
    chimera_seq=no_gap_sequence_from_alignment(''.join(ref_aln_list))
    indexed_ref_aln=[x for x in indexed_ref_aln if x!='-']
    chimera_index=indexed_ref_aln.index(aln_index)
    chi_AA=chimera_seq[chimera_index]
    return chimera_index,chi_AA,chimera_seq
S1='QCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPGSASS'
print(correct_alignment_for_chimera_index('6vsb_MSA.aln','6VSB',1767,'Myotis',S1))

def get_sequence_from_pdb(pdb_file, chain_id):
    # Iterate over all records in the PDB file
    for record in SeqIO.parse(pdb_file, "pdb-atom"):
        # Check if the record corresponds to the desired chain
        if record.annotations["chain"] == chain_id:
            # Return the sequence
            return str(record.seq)
    # Extract the sequence from the chain


def compare_aligned_contact(alignment_file, native_pdb, chimera_pdb, chain_id, label, alignment_index,
                            rank_difference_file=''):
    try:
        residue = get_residue_at_native_position(alignment_file, label, alignment_index)
        native_index = correct_alignment_for_residue_position(alignment_file, label, alignment_index)
    except:
        print('doesnt exist')
        return [label, alignment_index + 1, [], [], [], []]
    native_contacts = get_individual_intra_contacts(native_pdb, chain_id, native_index)
    native_seq = get_sequence_from_pdb(native_pdb, chain_id)
    if native_contacts:
        native_contacts_ids = [f'{native_seq[x]}{x + 1}' for x in native_contacts]
    else:
        native_contacts_ids = []
    native_contacts_ids += get_inter_protein_contacts(native_pdb, chain_id, True)[native_index]
    chimera_seq = get_sequence_from_pdb(chimera_pdb, chain_id)
    chi_index = correct_alignment_for_chimera_index(alignment_file, label, native_index, chimera_seq)
    residue_ids = f'{residue}{native_index + 1},{residue}{chi_index + 1}'
    chimera_contacts = get_individual_intra_contacts(chimera_pdb, chain_id, chi_index)
    if chimera_contacts:
        chimera_contacts_ids = [f'{chimera_seq[x]}{x + 1}' for x in chimera_contacts]
    else:
        chimera_contacts_ids = []
    chimera_contacts_ids += get_inter_protein_contacts(chimera_pdb, chain_id, True)[chi_index]
    if rank_difference_file:
        rank_difference = rank_difference_table_to_dict(rank_difference_file)[label]
        return [label, alignment_index + 1, residue_ids, rank_difference, native_contacts_ids, chimera_contacts_ids]
    else:
        return [label, alignment_index + 1, residue_ids, native_contacts_ids, chimera_contacts_ids]


def compare_aligned_contact_all_proteins(alignment_file, native_pdb_format, chimera_pdb_format, chain_id,
                                         alignment_index, label_list, contact_comparison_file, rank_difference_file=''):
    comparison_matrix = zeros((len(label_list), 5 + int(bool(rank_difference_file))), dtype=object)
    for index, label in enumerate(label_list):
        native = native_pdb_format.format(label)
        chi = chimera_pdb_format.format(label)
        comparison = compare_aligned_contact(alignment_file, native, chi, chain_id, label, alignment_index,
                                             rank_difference_file)
        print(comparison)
        comparison_matrix[index] = comparison
    convert_array_to_file(comparison_matrix, '/', contact_comparison_file)


aln = '/gpfs/gpfs0/scratch/jws6pq/Notebook/PDB/CoronavirusMSA.aln'
native_format = '/gpfs/gpfs0/scratch/jws6pq/Notebook/PDB/3mer{0}.pdb'
chi_format = '/gpfs/gpfs0/scratch/jws6pq/Notebook/PDB/3merSARS2w{0}S1.pdb'


# with open("/gpfs/gpfs0/scratch/jws6pq/Notebook/Overall/List_of_coronaviruses", 'r') as loc:
#     loc = [line.split()[-1] for line in loc.readlines()]
# compare_aligned_contact_all_proteins(aln,native_format,chi_format,'B',1267,loc,'/gpfs/gpfs0/scratch/jws6pq/Notebook/Immersion/CarbonB_1267_contacts.tsv','/gpfs/gpfs0/scratch/jws6pq/Notebook/Immersion/Rank_change_1268.tsv')
# def compare_aligned_contacts(alignment_file,reference_label,comparison_label,ref_pdb,com_pdb):
#     aligned_reference=ColorCoding.create_dictionary_from_alignment(alignment_file)[reference_label]
#     aligned_comparison = ColorCoding.create_dictionary_from_alignment(alignment_file)[comparison_label]
#     reference_indexing=get_alignment_indexing(aligned_reference)
#     comparison_indexing=get_alignment_indexing(aligned_comparison)
def ContactOverlap(alignment_file, comparison, reference='6VSB_B'):
    sequence_dictionary = create_dictionary_from_alignment(alignment_file)
    reference_sequence, comparison_sequence = sequence_dictionary[reference], sequence_dictionary[comparison]
    # CP is Comparison Protein and RP is Reference Protein
    cp_updated_contact_map = correct_contact_position_for_alignment(f'3mer{comparison}.pdb', 'B', comparison_sequence)
    print(comparison_sequence)
    print(cp_updated_contact_map)
    rp_updated_contact_map = correct_contact_position_for_alignment(f'{reference}.pdb', "B", reference_sequence)
    rp_contact_map, cp_contact_map = [], []
    j = 0
    residue_comparison = [blosum[where(blosum_62_matrix()[:, 0] == res_ref)[0][0] - 1][
                              where(blosum_62_matrix()[0, :] == res_com)[0][0] - 1]
                          for res_ref, res_com in zip(reference_sequence, comparison_sequence)]
    for residue in reference_sequence:
        if residue.isalpha():
            rp_contact_map.append(rp_updated_contact_map[j])
            j += 1
        else:
            rp_contact_map.append([])
    j = 0
    for residue in comparison_sequence:
        if residue.isalpha():
            cp_contact_map.append(cp_updated_contact_map[j])
            j += 1
        else:
            cp_contact_map.append([])
    fraction_conserved = []
    # print(list(zip(rp_contact_map,cp_contact_map)))

    for x, y in zip(rp_contact_map, cp_contact_map):
        if x == [] and y == []:
            fraction_conserved.append(1); continue
        elif x == [] or y == []:
            fraction_conserved.append(0); continue
        x_set = set();
        y_set = set()
        [x_set.update(num for num in range(contact - 6, contact + 7)) for contact in x]
        [y_set.update(num for num in range(contact - 6, contact + 7)) for contact in y]
        fraction_conserved.append((len([contact for contact in y if contact in x_set]) + len(
            [contact for contact in x if contact in y_set])) / (len(x) + len(y)))
    print(fraction_conserved)
