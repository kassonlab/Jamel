
def GetResidueContactPairs(PDBnickname,PDBFilename,DistanceCutoff):

    import Bio.PDB
    import numpy as np
    import sys
    pdb_code = PDBnickname
    pdb_filename = PDBFilename
    #Copied from https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/protein_contact_map/
    def calc_residue_dist(residue_one, residue_two) :
        """Returns the C-alpha distance between two residues"""
        diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord
        return np.sqrt(np.sum(diff_vector * diff_vector))

    # Copied from https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/protein_contact_map/
    def calc_dist_matrix(chain_one, chain_two) :
        """Returns a matrix of C-alpha distances between two chains"""
        answer = np.zeros((len(chain_one), len(chain_two)), np.float)
        for row, residue_one in enumerate(chain_one) :
            for col, residue_two in enumerate(chain_two) :
                if (calc_residue_dist(residue_one, residue_two))<=DistanceCutoff and (col-row)>=6:
                    answer[row, col] = 1
                else:
                    answer[row, col] = 0
        return answer

    structure = Bio.PDB.PDBParser().get_structure(pdb_code, pdb_filename)
    model = structure[0]
    np.set_printoptions(threshold=sys.maxsize)
    dist_matrix = calc_dist_matrix(model["B"], model["B"])
    #print(dist_matrix)
    X_axis=list(np.where(dist_matrix==1)[0])
    Y_axis=list(np.where(dist_matrix==1)[1])
    ListofContactPairs=[]
    i=0
    #print(zip(X_axis,Y_axis),1)
    for x,y in zip(X_axis,Y_axis):
        ListofContactPairs.append('('+str(x+1)+','+str(y+1)+')')
        i+=1
    print(ListofContactPairs)
GetResidueContactPairs('6vsb', '6VSB_B.pdb',8)
