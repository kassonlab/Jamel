import numpy as np
ProteinList=[line.split()[-1] for line in open('List','r').readlines()]
Domain='S1'
def FoldXStability(Protein,ComparisonScore=199.12):
    import os
    os.system('/sfs/lustre/bahamut/scratch/jws6pq/FoldX/foldx_20221231 --command=Stability --clean-mode=1 --pdb=SARS2w'+Protein+Domain+'.pdb --pdb-dir=/sfs/lustre/bahamut/scratch/jws6pq/Notebook/PDB --output-dir=/sfs/lustre/bahamut/scratch/jws6pq/Notebook/FoldXResults/')
    FoldxScore=float(open('/sfs/lustre/bahamut/scratch/jws6pq/Notebook/FoldXResults/SARS2w'+Protein+Domain+'_0_ST.fxout','r').readlines()[0].split()[1])
    FoldXDifference=-((ComparisonScore-FoldxScore)/ComparisonScore)*100
    return FoldXDifference
FoldXScore=list(map(FoldXStability,ProteinList))
def PieceWise(Protein,Domain='RBD',ComparisonProtein='alphasars2.pdb'):
    import RBDFinder
    cmd.load(ComparisonProtein,object='CP')
    cmd.remove('organic')
    proteinpdb = Protein + '.pdb'
    cmd.load(proteinpdb)
    # ComparisonProtein=CP
    cmd.select('CPNonDomain', selection='CP and not resi 1-233')
    cmd.select('CPDomain', selection='CP and resi -424')
    SpliceBoundaries=RBDFinder.AlignmentFinder(Protein+'onSARS2.aln', Protein)
    SpliceBoundary1=str(SpliceBoundaries[0])
    SpliceBoundary2=str(SpliceBoundaries[1])
    NTDindex = Protein+' and resi -'+SpliceBoundary1
    RBDindex = Protein+' and resi '+ SpliceBoundary1+'-'+SpliceBoundary2
    NTDName = Protein +'NTD'
    RBDName = Protein +'RBD'
    cmd.select(NTDName, selection=NTDindex)
    NTDRMSD=cmd.super('SARS2NTD', NTDName)[0]
    cmd.select(RBDName, selection=RBDindex)
    RBDRMSD = cmd.super('SARS2RBD', RBDName)[0]
    OverallRMSD = cmd.super('SARS2', Protein)[0]
    cmd.delete('all')
    return RBDRMSD, OverallRMSD
RMSD=list(map(PieceWise,ProteinList))
DomainRMSD=[x[0] for x in RMSD]
OverallRMSD=[x[1] for x in RMSD]
import Perresiduecomparison
import RBDFinder
AlignmentFileNames=[x+'onSARS2.aln' for x in ProteinList]
SpliceBoundaries=list(map(RBDFinder.AlignmentFinder,AlignmentFileNames,ProteinList))
ChimeraBoundary1=[x[0] for x in SpliceBoundaries]
ChimeraBoundary2=[x[1] for x in SpliceBoundaries]
SARS2Boundary1=[234 for x in ProteinList]
PercentDifference=list(map(Perresiduecomparison.ResidueConfidenceComparison,ProteinList,SARS2Boundary1,ChimeraBoundary1,ChimeraBoundary2))
OverallPercentDifference=list(map(Perresiduecomparison.OverallConfidenceComparison,ProteinList))
def SequenceSimilarity(Protein):
    EmbossScore=open(Protein+'.emboss','r').readlines()[27].split()[2]
    return EmbossScore

EmbossScore=list(map(SequenceSimilarity,ProteinList))

DataChart=np.empty((len(ProteinList)+1,7),dtype=object)
DataChart[0,0]='Protein'
DataChart[0,1]=Domain+'RMSD'
DataChart[0,2]='OverallRMSD'
DataChart[0,3]='Emboss'
DataChart[0,4]='DomainConfidence'
DataChart[0,5]='OverallConfidence'
DataChart[0,6]='FoldX'
DataChart[1:,0]=ProteinList
DataChart[1:,1]=RBDRMSD
DataChart[1:,2]=OverallRMSD
DataChart[1:,3]=EmbossScore
DataChart[1:,4]=PercentDifference
DataChart[1:,5]=OverallPercentDifference
DataChart[1:,6]=FoldXScore
np.savetxt('/sfs/lustre/bahamut/scratch/jws6pq/CMfiles/'+Domain+'ChimeraAnalysis.tsv', DataChart, fmt="%s,%s,%s,%s,%s,%s,%s", delimiter="")