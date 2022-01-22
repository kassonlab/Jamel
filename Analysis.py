import numpy as np

ProteinList=[line.split()[-1] for line in open('List','r').readlines()]
def PieceWise(Protein):
    import RBDFinder
    cmd.load("SARS2.pdb")
    cmd.remove('organic')
    proteinpdb = Protein + '.pdb'
    cmd.load(proteinpdb)
    cmd.select('SARS2NTD', selection='SARS2 and resi 1-233')
    cmd.select('SARS2RBD', selection='SARS2 and resi 234-432')
    cmd.select('SARS2Stalk', selection='SARS2 and resi 433-973')
    SpliceBoundaries=RBDFinder.AlignmentFinder(Protein+'onSARS2.aln', Protein)
    SpliceBoundary1=str(SpliceBoundaries[0])
    SpliceBoundary2=str(SpliceBoundaries[1])
    NTDindex = Protein+' and resi -'+SpliceBoundary1
    RBDindex = Protein+' and resi '+ SpliceBoundary1+'-'+SpliceBoundary2
    Stalkindex =Protein+ ' and resi '+ SpliceBoundary2+'-'
    NTDName = Protein +'NTD'
    RBDName = Protein +'RBD'
    StalkName = Protein +'Stalk'
    cmd.select(NTDName, selection=NTDindex)
    NTDRMSD=cmd.super('SARS2NTD', NTDName)[0]
    cmd.select(RBDName, selection=RBDindex)
    RBDRMSD = cmd.super('SARS2RBD', RBDName)[0]
    cmd.select(StalkName, selection=Stalkindex)
    StalkRMSD = cmd.super('SARS2Stalk', StalkName)[0]
    OverallRMSD = cmd.super('SARS2', Protein)[0]
    cmd.delete('all')
    return RBDRMSD, OverallRMSD
RMSD=list(map(PieceWise,ProteinList))
RBDRMSD=[x[0] for x in RMSD]
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
DataChart=np.empty((len(ProteinList)+1,6),dtype=object)
DataChart[0,0]='Protein'
DataChart[0,1]='RBDRMSD'
DataChart[0,2]='OverallRMSD'
DataChart[0,3]='Emboss'
DataChart[0,4]='RBDConfidence'
DataChart[0,5]='OverallConfidence'
DataChart[1:,0]=ProteinList
DataChart[1:,1]=RBDRMSD
DataChart[1:,2]=OverallRMSD
DataChart[1:,3]=EmbossScore
DataChart[1:,4]=PercentDifference
DataChart[1:,5]=OverallPercentDifference
np.savetxt('/sfs/lustre/bahamut/scratch/jws6pq/CMfiles/OverallChimeraAnalysis.tsv', DataChart, fmt="%s,%s,%s,%s,%s,%s", delimiter="")