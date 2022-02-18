def FoldXStability(Protein,Domain,ComparisonScore=199.12):
    import os
    os.system('/sfs/lustre/bahamut/scratch/jws6pq/FoldX/foldx_20221231 --command=Stability --clean-mode=1 --pdb=SARS2w'+Protein+Domain+'.pdb --pdb-dir=/sfs/lustre/bahamut/scratch/jws6pq/Notebook/PDB --output-dir=/sfs/lustre/bahamut/scratch/jws6pq/Notebook/FoldXResults/')
    FoldxScore=float(open('/sfs/lustre/bahamut/scratch/jws6pq/Notebook/FoldXResults/SARS2w'+Protein+Domain+'_0_ST.fxout','r').readlines()[0].split()[1])
    FoldXDifference=-((ComparisonScore-FoldxScore)/ComparisonScore)*100
    return FoldXDifference
def PieceWiseRMSD(Protein,CPSplice1,CPSplice2,SpliceBoundary1,SpliceBoundary2,ComparisonProtein='6VSB_B.pdb'):
    import RBDFinder
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
def SpliceConfidenceComparison(Protein,ChimeraSplice1,ChimeraSplice2,SARS2Splice1,Domain):
    ProteinScore=list(map(float,open(Protein+'.plddt', 'r').readlines()))
    ChimeraScore=list(map(float,open('SARS2w'+Protein+Domain+'.plddt', 'r').readlines()))
    SpliceLength=ChimeraSplice2-ChimeraSplice1
    OriginalAverageScore=sum(ProteinScore[(ChimeraSplice1):(ChimeraSplice2)])/(SpliceLength)
    ChimeraAverageScore=sum(ChimeraScore[(int(SARS2Splice1)-1):(SpliceLength+int(SARS2Splice1))])/(SpliceLength)
    ScorePercentDifference=(OriginalAverageScore-ChimeraAverageScore)/OriginalAverageScore
    return ScorePercentDifference*100
def OverallConfidenceComparison(Protein,Domain):
    SARS2Score = list(map(float, open('SARS2.plddt', 'r').readlines()))
    ChimeraScore = list(map(float, open('SARS2w' + Protein + Domain+'.plddt', 'r').readlines()))
    SARS2AverageScore = sum(SARS2Score) / len(SARS2Score)
    ChimeraAverageScore = sum(ChimeraScore) / len(ChimeraScore)
    ScorePercentDifference = (SARS2AverageScore - ChimeraAverageScore) / SARS2AverageScore
    return ScorePercentDifference * 100
print(SequenceSimilarity('ShorebirdDelta','RBD'))