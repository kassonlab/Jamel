from MDAnalysis.analysis import distances as MDdist
from MDAnalysis import Universe # test trajectory
from numpy import any as npany
from numpy import save,load,sqrt, array, where,  zeros
from os.path import exists
import json
from Bio import PDB, SeqIO
from Analysis import convert_array_to_file
# TODO allow looking at contacts for specific position for intra
# TODO combine contacts across all homomers
from scipy.spatial import KDTree
from AccessiontoAlignment import get_alignment_indexing,get_alignment_indexing_w_dashes, alignment_finder,create_dictionary_from_alignment,no_gap_sequence_from_alignment

#New MDAnalysis stuff
# TODO for kdtree
def intra_residue_dist_matrix_md_analysis(pdb_file, chain_id,residue_separation=6, npy_matrix_file='', distance_cutoff=6):
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
        save(npy_matrix_file,residue_dist_matrix)
    return residue_dist_matrix
def individual_intra_residue(npy_matrix_file,residue_index):
    residue_row=load(npy_matrix_file)[residue_index,:]

    return [ind for ind,x in enumerate(residue_row) if x==1]
def check_for_contact(atom_group_1,atom_group_2,universe,distance_cutoff=6):
    dist_arr = MDdist.distance_array(atom_group_1,  # reference
                                     atom_group_2,  # configuration
                                     box=universe.dimensions)
    return npany(dist_arr < distance_cutoff)
def segment_to_str(segment_obj):
    return str(segment_obj).split()[-1].replace('>', "")


def inter_residue_contact_list_md_analysis(pdb_file, chain_id, saved_contact_json='', distance_cutoff=6):
    """Takes a distance matrix of a protein and converts it to a binary matrix that uses 1 to signify a residue
    contact and 0 for no contact or Returns a matrix of C-alpha distances between residues in a protein chain or b"""
    if saved_contact_json and exists(saved_contact_json):
        with open(saved_contact_json, 'r') as file:
            loaded_list = json.load(file)
        inter_contacts=loaded_list
        return inter_contacts
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
                if check_for_contact(res,foreign_res,protein_coords,distance_cutoff):
                    inter_contacts[index].append(f'{foreign_chain}:{foreign_index}')
    if saved_contact_json and not exists(saved_contact_json):
        with open(saved_contact_json, 'w') as file:
            json.dump(inter_contacts, file)
    return inter_contacts

def get_intra_residue_contact_list(pdb_filename, chain_identifier, npy_matrix_file):
    """Takes a pdb structure and return a nested list where each residue index in the list has a list with the python
    indexes for the residues it's in contact with.
    Residue indexes that have no contacts will have an empty list"""
    if exists(npy_matrix_file+'.npy'):
        dist_matrix =load(npy_matrix_file+'.npy')
    else:
        dist_matrix = intra_residue_dist_matrix_md_analysis(pdb_filename,chain_identifier,npy_matrix_file=npy_matrix_file)
    x_axis,y_axis=list(where(dist_matrix==1)[0]),list(where(dist_matrix==1)[1])
    list_of_contact_pairs=[[] for rows in dist_matrix]
    for x, y in zip(x_axis, y_axis):
        list_of_contact_pairs[x].append(y)
    return list_of_contact_pairs




