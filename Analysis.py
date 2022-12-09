import os


def RunFoldX(Protein,Domain,ComparisonScore=827.13):
    import os
    os.system('/scratch/jws6pq/FoldX/foldx_20221231 --command=Stability --clean-mode=1 --pdb=SARS2w'+Protein+Domain+'.pdb --pdb-dir=/scratch/jws6pq/Notebook/PDB --output-dir=/scratch/jws6pq/Notebook/FoldXResults/')
def FoldXStability(Protein,Domain,ComparisonScore=827.13):
    FoldxScore=float(open('/scratch/jws6pq/Notebook/FoldXResults/SARS2w'+Protein+Domain+'_0_ST.fxout','r').readlines()[0].split()[1])
    FoldXDifference=-((ComparisonScore-FoldxScore)/ComparisonScore)*100
    return FoldXDifference
def PieceWiseRMSD(Protein,CPSplice1,CPSplice2,SpliceBoundary1,SpliceBoundary2,ComparisonProtein='6VSB_B.pdb'):
    from pymol import cmd
    cmd.load(ComparisonProtein,object='CP')
    cmd.remove('organic')
    proteinpdb = Protein + '.pdb'
    cmd.load(proteinpdb)
    # ComparisonProtein=CP
    cmd.select('CPNonDomain', selection='CP and not resi '+str(CPSplice1)+'-'+str(CPSplice2))
    cmd.select('CPDomain', selection='CP and resi '+str(CPSplice1)+'-'+str(CPSplice2))
    Domainindex = Protein+' and resi '+str(SpliceBoundary1)+'-'+str(SpliceBoundary2)
    NonDomainindex = Protein+' and not resi '+ str(SpliceBoundary1)+'-'+str(SpliceBoundary2)
    DomainName = Protein +'Domain'
    NonDomainName = Protein +'NonDomain'
    cmd.select(DomainName, selection=Domainindex)
    cmd.select(NonDomainName, selection=NonDomainindex)
    DomainRMSD=cmd.align('CPDomain', DomainName)[0]
    cmd.delete('all')
    return DomainRMSD
def OverallRMSD(Protein,Comparison='3merSARS2.pdb'):
    from pymol import cmd
    proteinpdb = Protein
    cmd.load(Comparison, object='CP')
    cmd.load(proteinpdb,object='Protein')
    cmd.remove('organic')
    cmd.load(proteinpdb)
    RMSD = cmd.align('CP', 'Protein')[0]
    cmd.delete('all')
    return RMSD

def SequenceSimilarity(Protein,Domain):
    EmbossScore=open(Protein+Domain+'.emboss','r').readlines()[25].split()[-1]
    return EmbossScore.replace('(','').replace(')','').replace('%','')
def OverallConfidence(plddtfile):
    plddt= list(map(float, open(plddtfile, 'r').readlines()))
    Averageplddt=sum(plddt)/len(plddt)
    return Averageplddt


def ConfidenceComparison(Protein,ChimeraSplice1,ChimeraSplice2,SARS2Splice1,SARS2Splice2,Domain):
    SARS2Score = list(map(float, open('SARS2.plddt', 'r').readlines()))
    ProteinScore=list(map(float, open(Protein+'.plddt', 'r').readlines()))
    ChimeraScore=list(map(float,open('SARS2w'+Protein+Domain+'.plddt', 'r').readlines()))
    j=0
    k=0
    SpliceLength=ChimeraSplice2-ChimeraSplice1
    ScoreDifference=0
    ChimeraLength=len(ChimeraScore)
    for ChimeraResidue in range(ChimeraLength):
        if SARS2Splice1<=ChimeraResidue<=(SARS2Splice1+SpliceLength):
            ScoreDifference+=(ChimeraScore[ChimeraResidue]-ProteinScore[ChimeraSplice1+j])/ProteinScore[ChimeraSplice1+j]
            j+=1
        elif 0<=ChimeraResidue<SARS2Splice1:
            #This statement is calculating the difference for the first part of SARS where there was no cleaving
            ScoreDifference+=(ChimeraScore[ChimeraResidue]-SARS2Score[ChimeraResidue])/SARS2Score[ChimeraResidue]
        else:
            ScoreDifference+=(ChimeraScore[ChimeraResidue]-SARS2Score[SARS2Splice2+1+k])/SARS2Score[SARS2Splice2+1+k]
            k+=1
    AveragePercentScoreDifference=ScoreDifference/ChimeraLength*100

    return AveragePercentScoreDifference
def MultimerConfidenceComparison(nativeplddt,Chimeraplddt,Chimeraboundarytuple,NativeBoundaryTuple):
    NativeProteinScore = list(map(float, open(nativeplddt, 'r').readlines()))
    ChimeraScore=list(map(float,open(Chimeraplddt, 'r').readlines()))
    SpliceLength=len(ChimeraScore[Chimeraboundarytuple[0]:Chimeraboundarytuple[1]])
    Relativedifference=0
    Nativerange=[ind+NativeBoundaryTuple[0] for ind,x in enumerate(NativeProteinScore[NativeBoundaryTuple[0]:NativeBoundaryTuple[1]])]
    Chimerarange = [ind+Chimeraboundarytuple[0] for ind,x in enumerate(ChimeraScore[Chimeraboundarytuple[0]:Chimeraboundarytuple[1]])]
    for x,y in zip(Chimerarange,Nativerange):
        Relativedifference+=(ChimeraScore[x]-NativeProteinScore[y])/NativeProteinScore[y]*100
    Relativedifference=Relativedifference

    return Relativedifference,SpliceLength

def AveragingMultimerPLDDT(Plddtfilename,Subunits=3):
    MultimerPlddt=list(map(float, open(Plddtfilename, 'r').readlines()))
    Monomerlength=int(len(MultimerPlddt)/int(Subunits))
    Newplddtfile = open('Avg'+Plddtfilename, 'w')
    ResidueIndex=0
    while ResidueIndex in range(Monomerlength):
        SubunitIndex = 0
        AveragedPlddt=0
        #You can step through list
        while SubunitIndex in range(Subunits):
            AveragedPlddt+=MultimerPlddt[ResidueIndex+(Monomerlength*SubunitIndex)]
            SubunitIndex += 1
        AveragedPlddt=AveragedPlddt/Subunits
        Newplddtfile.write(str(AveragedPlddt) + '\n')
        ResidueIndex+=1
    Newplddtfile.close()
    return Newplddtfile.name
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
    EmbossScore = open('/gpfs/gpfs0/scratch/jws6pq/Notebook/Emboss/Full' + comparison + '.emboss', 'r').readlines()[25].split()[-1]
    return comparison,TotalContacts,EmbossScore.replace('(','').replace(')','').replace('%','')
#Do i consider all the times where there are residues beyond SARS???????
def FaultScan(proteinpdb):
    return 1 if OverallRMSD(proteinpdb)>35.5 else 0
def RMSF(Protein,Timestepinps):
    from os import system as sys
    sys('echo 1 | gmx_mpi trjconv -f /gpfs/gpfs0/scratch/jws6pq/Gromacs/'+Protein+'.xtc -s /gpfs/gpfs0/scratch/jws6pq/Gromacs/3merSARS2wBatHKU4S1_production_1.tpr -dt '+Timestepinps+' -o /scratch/jws6pq/Gromacs/'+Protein+'RMSF.xtc')
    sys('echo 1 | gmx_mpi rmsf -dt '+Timestepinps+' -res -f /gpfs/gpfs0/scratch/jws6pq/Gromacs/'+Protein+'RMSF.xtc -s /gpfs/gpfs0/scratch/jws6pq/Gromacs/'+Protein+'_production_1.tpr -o '+Protein+'rmsf.xvg')