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
    NTDRMSD=cmd.align('SARS2NTD', NTDName)[0]
    cmd.select(RBDName, selection=RBDindex)
    RBDRMSD = cmd.align('SARS2RBD', RBDName)[0]
    cmd.select(StalkName, selection=Stalkindex)
    StalkRMSD = cmd.align('SARS2Stalk', StalkName)[0]
    OverallRMSD = cmd.align('SARS2', Protein)[0]
    return RBDRMSD
import numpy as np
RMSD=list(map(PieceWise,ProteinList))
print(RMSD)
import Perresiduecomparison
import RBDFinder
AlignmentFileNames=[x+'onSARS2.aln' for x in ProteinList]
SpliceBoundaries=list(map(RBDFinder.AlignmentFinder,AlignmentFileNames,ProteinList))
ChimeraBoundary1=[x[0] for x in SpliceBoundaries]
ChimeraBoundary2=[x[1] for x in SpliceBoundaries]
SARS2Boundary1=[234 for x in ProteinList]
SARS2Boundary2=[432 for x in ProteinList]
PercentDifference=list(map(Perresiduecomparison.ConfidenceComparison,ProteinList,SARS2Boundary1,SARS2Boundary2,ChimeraBoundary1,ChimeraBoundary2))
print(PercentDifference)
def SequenceSimilarity(Protein):
    EmbossScore=open(Protein+'.emboss','r').readlines()[27].split()[2]
    return EmbossScore
EmbossScore=list(map(SequenceSimilarity,ProteinList))
print(EmbossScore)
