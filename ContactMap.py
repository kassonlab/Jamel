def calc_dist_matrix(chain_one, chain_two,DistanceCutoff):
    from numpy import array
    # """Returns a matrix of C-alpha distances between two chains"""
    answer=array([1 if (calc_residue_dist(residue_one, residue_two)) <= DistanceCutoff and (col-row) >= 6 else 0 for row, residue_one in enumerate(chain_one) for col, residue_two in enumerate(chain_two)]).reshape((len(chain_one), len(chain_two)))
    return array(answer)
#Copied from https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/protein_contact_map/
def calc_residue_dist(residue_one, residue_two) :
    from numpy import sqrt
        # """Returns the C-alpha distance between two residues"""
    diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord
    return sqrt(sum(diff_vector * diff_vector))
def GetResidueContactPairs(PDBnickname,PDBFilename,DistanceCutoff):
    from Bio.PDB import PDBParser
    from numpy import where
    pdb_code,pdb_filename = PDBnickname,PDBFilename
    structure = PDBParser().get_structure(pdb_code, pdb_filename)
    model = structure[0]
    dist_matrix = calc_dist_matrix(model["B"], model["B"],DistanceCutoff)
    X_axis,Y_axis=list(where(dist_matrix==1)[0]),list(where(dist_matrix==1)[1])
    ListofContactPairs=[[] for x in dist_matrix]
    for x, y in zip(X_axis, Y_axis):
        ListofContactPairs[x].append(y)
    return ListofContactPairs
def CorrectResiduePositionforAlignment(Protein,alignment):
    ContactMap = GetResidueContactPairs(Protein, '3mer' + Protein + '.pdb', 7)
    SeqIndexing = [ind for ind, x in enumerate(alignment) if x != '-']
    ResiduePositionDictionary = {indx: indy for indx, indy in enumerate(SeqIndexing)}
    UpdatedContactMap = [[ResiduePositionDictionary[y] for y in ContactMap[ind]] for ind, x in enumerate(ContactMap)]
    return UpdatedContactMap
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
    # EmbossScore = open('/gpfs/gpfs0/scratch/jws6pq/Notebook/Emboss/Full' + comparison + '.emboss', 'r').readlines()[25].split()[-1]
    # return comparison,Overlap,EmbossScore.replace('(','').replace(')','').replace('%','')
#Do i consider all the times where there are residues beyond SARS???????
ContactOverlap('SARS2wEverythingstable.aln','WhiteHKU16')



