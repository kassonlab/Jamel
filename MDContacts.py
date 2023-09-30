from MDAnalysis.analysis import distances as MDdist
from MDAnalysis import Universe # test trajectory
from numpy import any as npany
from numpy import zeros,save,load,where
from ContactMap import get_intra_residue_contact_pairs
#New MDAnalysis stuff
# target_chain = "A"  # Replace "A" with the chain ID you want to select
#
# # Select atoms from the specified chain
# atoms_from_target_chain = u.select_atoms(f"chain {target_chain}")
# TODO for kdtree
def intra_residue_dist_matrix_md_analysis(pdb_file, chain_id,residue_separation=6, npy_matrix_file='',make_it_binary='Yes', distance_cutoff=6):
    """Takes a distance matrix of a protein and converts it to a binary matrix that uses 1 to signify a residue
    contact and 0 for no contact or Returns a matrix of C-alpha distances between residues in a protein chain or b"""
    protein_coords = Universe(pdb_file, pdb_file)
    # Turns residues into atom groups and excludes non heavy atoms (hydrogens)
    atom_groups = [res.atoms.select_atoms(f'not type H and chainID {chain_id}') for res in protein_coords.residues if len(res.atoms.select_atoms(f'not type H and chainID {chain_id}'))>0]
    residue_dist_matrix=zeros((len(atom_groups),len(atom_groups)))
    for y_index,res in enumerate(atom_groups):
        for x_index,res in enumerate(atom_groups):
            # WHat residue separation should I use
            if x_index not in range(y_index-residue_separation,y_index+residue_separation):
                residue_dist_matrix[y_index,x_index]=int(check_for_contact(atom_groups[y_index],atom_groups[x_index],protein_coords,distance_cutoff))
    if npy_matrix_file:
        save(npy_matrix_file, residue_dist_matrix)
    return residue_dist_matrix

def check_for_contact(atom_group_1,atom_group_2,universe,distance_cutoff=6):
    dist_arr = MDdist.distance_array(atom_group_1,  # reference
                                     atom_group_2,  # configuration
                                     box=universe.dimensions)
    return npany(dist_arr < distance_cutoff)
def segment_to_str(segment_obj):
    return str(segment_obj).split()[-1].replace('>', "")


def inter_residue_dist_matrix_md_analysis(pdb_file, chain_id,residue_separation=6, npy_matrix_file='',make_it_binary='Yes', distance_cutoff=6):
    """Takes a distance matrix of a protein and converts it to a binary matrix that uses 1 to signify a residue
    contact and 0 for no contact or Returns a matrix of C-alpha distances between residues in a protein chain or b"""
    protein_coords = Universe(pdb_file, pdb_file)
    # Turns residues into atom groups and excludes non heavy atoms (hydrogens)
    chain_dict={segment_to_str(chain_id):[] for chain_id in protein_coords.segments}
    for res in protein_coords.residues:
        chain_dict[segment_to_str(res.segment)].append(res.atoms.select_atoms(f'not type H'))
    residues_of_interest=chain_dict[chain_id]
    del chain_dict[chain_id]
    inter_contacts=[[] for res in residues_of_interest]
    for index,res in enumerate(residues_of_interest):
        for foreign_chain,res_list in chain_dict.items():
            for foreign_index,foreign_res in enumerate(res_list):
                if check_for_contact(res,foreign_res,protein_coords):
                    inter_contacts[index].append(f'{foreign_chain}:{foreign_index}')
    return inter_contacts

def get_residue_contact_pairs(pdb_filename,chain_identifier,npy_matrix_file=''):
    """Takes a pdb structure and return a nested list where each residue index in the list has a list with the python
    indexes for the residues it's in contact with.
    Residue indexes that have no contacts will have an empty list"""
    if npy_matrix_file:
        dist_matrix =load(npy_matrix_file)
    else:
        dist_matrix = intra_residue_dist_matrix_md_analysis(pdb_filename,chain_identifier,npy_matrix_file=npy_matrix_file)
    x_axis,y_axis=list(where(dist_matrix==1)[0]),list(where(dist_matrix==1)[1])
    list_of_contact_pairs=[[] for rows in dist_matrix]
    for x, y in zip(x_axis, y_axis):
        list_of_contact_pairs[x].append(y)
    return list_of_contact_pairs
# turn matrix into a file npy file prolly
print(inter_residue_dist_matrix_md_analysis('3merSARS2.pdb','A'))
