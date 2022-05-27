
    import Bio.PDB
    import numpy as np

    pdb_code = "6VSB"
    pdb_filename = "6VSB_B.pdb"

    def calc_residue_dist(residue_one, residue_two) :
        """Returns the C-alpha distance between two residues"""
        diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord
        return np.sqrt(np.sum(diff_vector * diff_vector))

    def calc_dist_matrix(chain_one, chain_two) :
        """Returns a matrix of C-alpha distances between two chains"""
        answer = np.zeros((len(chain_one), len(chain_two)), np.float)

        for row, residue_one in enumerate(chain_one) :
            for col, residue_two in enumerate(chain_two) :
                answer[row, col] = calc_residue_dist(residue_one, residue_two)
        return answer

    structure = Bio.PDB.PDBParser().get_structure(pdb_code, pdb_filename)
    model = structure[0]

    dist_matrix = calc_dist_matrix(model["B"], model["B"])
    # file = open('ContactMap', 'w')
    # file.write(dist_matrix)
    # file.close()
    X_axis=list(np.where(dist_matrix<=7)[0])
    Y_axis=list(np.where(dist_matrix<=7)[1])
    print(dist_matrix[7,149])
    print(X_axis)
    print(Y_axis)
    ListofContactPairs=[]
    i=0
    for x,y in zip(X_axis,Y_axis):
        if y-x>5:
            ListofContactPairs.append('('+str(x+1)+','+str(y+1)+')')
            i+=1

print(ListofContactPairs[1296],ListofContactPairs)
# np.savetxt('Contct',dist_matrix,fmt="%s",delimiter=" ")