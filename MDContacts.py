from MDAnalysis.analysis import distances as MDdist
from MDAnalysis import Universe # test trajectory
from numpy import any as npany
from numpy import zeros
#New MDAnalysis stuff
# target_chain = "A"  # Replace "A" with the chain ID you want to select
#
# # Select atoms from the specified chain
# atoms_from_target_chain = u.select_atoms(f"chain {target_chain}")
# TODO for kdtree
def intra_residue_dist_matrix_md_analysis(pdb_file, chain_id,residue_separation=6, make_it_binary='Yes', distance_cutoff=6):
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
                dist_arr = MDdist.distance_array(atom_groups[y_index].positions,  # reference
                                                 atom_groups[x_index].positions,  # configuration
                                                 box=protein_coords.dimensions)
                if npany(dist_arr < distance_cutoff):
                    residue_dist_matrix[y_index,x_index]=1
    return residue_dist_matrix

# turn matrix into a file npy file prolly
print(intra_residue_dist_matrix_md_analysis('3merSARS2.pdb','B')[0])