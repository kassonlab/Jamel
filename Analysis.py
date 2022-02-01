import numpy as np
import Perresiduecomparison
import RBDFinder
def FoldXStability(Protein,Domain='RBD',ComparisonScore=199.12):
    import os
    os.system('/sfs/lustre/bahamut/scratch/jws6pq/FoldX/foldx_20221231 --command=Stability --clean-mode=1 --pdb=SARS2w'+Protein+Domain+'.pdb --pdb-dir=/sfs/lustre/bahamut/scratch/jws6pq/Notebook/PDB --output-dir=/sfs/lustre/bahamut/scratch/jws6pq/Notebook/FoldXResults/')
    FoldxScore=float(open('/sfs/lustre/bahamut/scratch/jws6pq/Notebook/FoldXResults/SARS2w'+Protein+Domain+'_0_ST.fxout','r').readlines()[0].split()[1])
    FoldXDifference=-((ComparisonScore-FoldxScore)/ComparisonScore)*100
    return FoldXDifference
def PieceWise(Protein,CPSplice1='224',CPSplice2='424',Domain='RBD',ComparisonProtein='alphasars2.pdb'):
    import RBDFinder
    cmd.load(ComparisonProtein,object='CP')
    cmd.remove('organic')
    proteinpdb = Protein + '.pdb'
    cmd.load(proteinpdb)
    # ComparisonProtein=CP
    cmd.select('CPNonDomain', selection='CP and not resi '+CPSplice1+'-'+CPSplice2)
    cmd.select('CPDomain', selection='CP and resi '+CPSplice1+'-'+CPSplice2)
    SpliceBoundaries=RBDFinder.AlignmentFinder(Protein+'onSARS2.aln', Protein)
    SpliceBoundary1=str(SpliceBoundaries[0])
    SpliceBoundary2=str(SpliceBoundaries[1])
    Domainindex = Protein+' and resi '+SpliceBoundary1+'-'+SpliceBoundary2
    NonDomainindex = Protein+' and not resi '+ SpliceBoundary1+'-'+SpliceBoundary2
    DomainName = Protein +'Domain'
    NonDomainName = Protein +'NonDomain'
    cmd.select(DomainName, selection=Domainindex)
    cmd.select(NonDomainName, selection=NonDomainindex)
    DomainRMSD=cmd.super('CPDomain', DomainName)[0]
    NonDomainRMSD = cmd.super('CPNonDomain', NonDomainName)[0]
    OverallRMSD = cmd.super('CP', Protein)[0]
    cmd.delete('all')
    return DomainRMSD, OverallRMSD
def SequenceSimilarity(Protein):
    EmbossScore=open(Protein+'.emboss','r').readlines()[27].split()[2]
    return EmbossScore
ProteinList=[line.split()[-1] for line in open('List','r').readlines()]
AlignmentFileNames=[x+'onSARS2.aln' for x in ProteinList]
SpliceBoundaries=list(map(RBDFinder.AlignmentFinder,AlignmentFileNames,ProteinList))
ChimeraBoundary1=[x[0] for x in SpliceBoundaries]
ChimeraBoundary2=[x[1] for x in SpliceBoundaries]
DomainSetting='S1'
RMSD=list(map(PieceWise,ProteinList,Domain=DomainSetting))
DomainRMSD=[x[0] for x in RMSD]
OverallRMSD=[x[1] for x in RMSD]
PercentDifference=list(map(Perresiduecomparison.ResidueConfidenceComparison,ProteinList,SARS2Boundary1,ChimeraBoundary1,ChimeraBoundary2))
OverallPercentDifference=list(map(Perresiduecomparison.OverallConfidenceComparison,ProteinList))
EmbossScore=list(map(SequenceSimilarity,ProteinList))
FoldXScore=list(map(FoldXStability,ProteinList,Domain=DomainSetting))
DataChart=np.empty((len(ProteinList)+1,7),dtype=object)
DataChart[0,0]='Protein'
DataChart[0,1]=DomainSetting+'RMSD'
DataChart[0,2]='OverallRMSD'
DataChart[0,3]='Emboss'
DataChart[0,4]='DomainConfidence'
DataChart[0,5]='OverallConfidence'
DataChart[0,6]='FoldX'
DataChart[1:,0]=ProteinList
DataChart[1:,1]=DomainRMSD
DataChart[1:,2]=OverallRMSD
DataChart[1:,3]=EmbossScore
DataChart[1:,4]=PercentDifference
DataChart[1:,5]=OverallPercentDifference
DataChart[1:,6]=FoldXScore
np.savetxt('/sfs/lustre/bahamut/scratch/jws6pq/CMfiles/'+DomainSetting+'ChimeraAnalysis.tsv', DataChart, fmt="%s,%s,%s,%s,%s,%s,%s", delimiter="")
