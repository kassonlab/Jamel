
def calc_dist_matrix(chain_one, chain_two,DistanceCutoff):
    import numpy as np
    # """Returns a matrix of C-alpha distances between two chains"""
    answer = np.zeros((len(chain_one), len(chain_two)), float)
    for row, residue_one in enumerate(chain_one):
        for col, residue_two in enumerate(chain_two):
            if (calc_residue_dist(residue_one, residue_two)) <= DistanceCutoff and (col - row) >= 6:
                answer[row, col] = 1
            else:
                answer[row, col] = 0
    return answer
#Copied from https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/protein_contact_map/
def calc_residue_dist(residue_one, residue_two) :
    import numpy as np
        # """Returns the C-alpha distance between two residues"""
    diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord
    return np.sqrt(np.sum(diff_vector * diff_vector))
def GetResidueContactPairs(PDBnickname,PDBFilename,DistanceCutoff):

    import Bio.PDB
    import numpy as np
    import sys
    pdb_code = PDBnickname
    pdb_filename = PDBFilename
    structure = Bio.PDB.PDBParser().get_structure(pdb_code, pdb_filename)
    model = structure[0]
    # np.set_printoptions(threshold=sys.maxsize)
    dist_matrix = calc_dist_matrix(model["B"], model["B"],DistanceCutoff)
    X_axis=list(np.where(dist_matrix==1)[0])
    Y_axis=list(np.where(dist_matrix==1)[1])
    ListofContactPairs=[[] for x in dist_matrix]
    i=0
    for x, y in zip(X_axis, Y_axis):
        ListofContactPairs[x].append(y)
    return ListofContactPairs
def ContactOverlap(Alignmentfile):
    # Alignment in FASTA format. Make sure your benchmark sequence is first
    Sequences = open(Alignmentfile, "r").read().split('>')
    SequenceDictionary={item.split('\n')[0]:item.split('\n')[1].strip() for item in Sequences if len(item)!=0}
    ScoreArray=np.empty((len(SequenceDictionary)-1,2), dtype=object)
    j=0
    for key, value in SequenceDictionary.items():
        ContactMap=GetResidueContactPairs(key,'3mer'+key+'.pdb',7)
        SeqIndexing = [ind for ind, x in enumerate(value) if x != '-']
        ResiduePositionDictionary={indx:indy for indx,indy in enumerate(SeqIndexing)}
        i=0
        UpdatedContactMap=[]
        for x in value:
            if x=='-':
                UpdatedContactMap.append([x])
            else:
                MapList = [x]
                for y in ContactMap[i]:
                    MapList+=[ResiduePositionDictionary[y]]
                i+=1
                UpdatedContactMap.append(MapList)
        if key=='6vsb_B':
            SARSContactMap=UpdatedContactMap
        else:
            Difference=0
            for x,y in zip(SARSContactMap,UpdatedContactMap):
                for w in x:
                    if w in y and w!='-':
                        Difference+=1
            ScoreArray[j,0]=key
            ScoreArray[j,1]=Difference
    np.savetxt('ContactScore.tsv', ScoreArray, fmt="%s,%s", delimiter="")
    return ScoreArray
#Do i consider all the times where there are residues beyond SARS???????


ContactOverlap('SARS2wEverythingstable.aln')
