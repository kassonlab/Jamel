def RunFoldX(Protein,Domain,ComparisonScore=827.13):
    import os
    os.system('/sfs/lustre/bahamut/scratch/jws6pq/FoldX/foldx_20221231 --command=Stability --clean-mode=1 --pdb=SARS2w'+Protein+Domain+'.pdb --pdb-dir=/sfs/lustre/bahamut/scratch/jws6pq/Notebook/PDB --output-dir=/sfs/lustre/bahamut/scratch/jws6pq/Notebook/FoldXResults/')
def FoldXStability(Protein,Domain,ComparisonScore=827.13):
    FoldxScore=float(open('/sfs/lustre/bahamut/scratch/jws6pq/Notebook/FoldXResults/SARS2w'+Protein+Domain+'_0_ST.fxout','r').readlines()[0].split()[1])
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
    DomainRMSD=cmd.super('CPDomain', DomainName)[0]
    NonDomainRMSD = cmd.super('CPNonDomain', NonDomainName)[0]
    OverallRMSD = cmd.super('CP', Protein)[0]
    cmd.delete('all')
    return DomainRMSD, OverallRMSD
def SequenceSimilarity(Protein,Domain):
    EmbossScore=open(Protein+Domain+'.emboss','r').readlines()[25].split()[-1]
    return EmbossScore.replace('(','').replace(')','').replace('%','')
def ConfidenceComparison(Protein,ChimeraSplice1,ChimeraSplice2,SARS2Splice1,SARS2Splice2,Domain):

    SARS2Score = list(map(float, open('SARS2.plddt', 'r').readlines()))
    ProteinScore=list(map(float, open(Protein+'.plddt', 'r').readlines()))
    ChimeraScore=list(map(float,open('SARS2w'+Protein+Domain+'.plddt', 'r').readlines()))
    ChimeraResidue=0
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
def MultimerConfidenceComparison(Protein,ChimeraSplice1,ChimeraSplice2,ComparisonProteinSplice1,ComparisonProteinSplice2,Domain,ComparisonProtein):
    ComparisonProteinScore = list(map(float, open('Avg'+ComparisonProtein+'.plddt', 'r').readlines()))
    ProteinScore=list(map(float, open('Avg'+Protein+'.plddt', 'r').readlines()))
    Protein=Protein[4:]
    ChimeraScore=list(map(float,open('Avg'+ComparisonProtein+'w'+Protein+Domain+'.plddt', 'r').readlines()))
    ChimeraResidue=0
    j=0
    k=0
    SpliceLength=ChimeraSplice2-ChimeraSplice1
    ScoreDifference=0
    ChimeraLength=len(ChimeraScore)
    for ChimeraResidue in range(ChimeraLength):
        if ComparisonProteinSplice1<=ChimeraResidue<=(ComparisonProteinSplice1+SpliceLength):
            ScoreDifference+=(ChimeraScore[ChimeraResidue]-ProteinScore[ChimeraSplice1+j])/ProteinScore[ChimeraSplice1+j]
            j+=1
        elif 0<=ChimeraResidue<ComparisonProteinSplice1:
            #This statement is calculating the difference for the first part of SARS where there was no cleaving
            ScoreDifference+=(ChimeraScore[ChimeraResidue]-ComparisonProteinScore[ChimeraResidue])/ComparisonProteinScore[ChimeraResidue]
        else:
            ScoreDifference+=(ChimeraScore[ChimeraResidue]-ComparisonProteinScore[ComparisonProteinSplice2+1+k])/ComparisonProteinScore[ComparisonProteinSplice2+1+k]

            k+=1

    AveragePercentScoreDifference=ScoreDifference/ChimeraLength*100
    return AveragePercentScoreDifference
#make sure the number of multimers is indicated at the front of the filename
def AveragingMultimerPLDDT(Plddtfilename):
    MultimerPlddt=list(map(float, open(Plddtfilename, 'r').readlines()))
    Subunits=int(Plddtfilename[0])
    Monomerlength=int(len(MultimerPlddt)/int(Subunits))
    Newplddtfile = open('Avg'+Plddtfilename, 'w')
    ResidueIndex=0
    while ResidueIndex in range(Monomerlength):
        SubunitIndex = 0
        AveragedPlddt=0
        while SubunitIndex in range(Subunits-1):
            AveragedPlddt+=MultimerPlddt[ResidueIndex+(Monomerlength*SubunitIndex)]
            AveragedPlddt=AveragedPlddt/Subunits
            SubunitIndex+=1
        Newplddtfile.write(str(AveragedPlddt) + '\n')
        ResidueIndex+=1
    Newplddtfile.close()
    return Newplddtfile.name
