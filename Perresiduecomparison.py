def ResidueConfidenceComparison(Protein,Boundary1,Boundary3,Boundary4):
    ProteinScore=list(map(float,open(Protein+'.plddt', 'r').readlines()))
    ChimeraScore=list(map(float,open('SARS2w'+Protein+'RBD.plddt', 'r').readlines()))
    SpliceLength=Boundary4-Boundary3
    OriginalAverageScore=sum(ProteinScore[(Boundary3-1):(Boundary4)])/(SpliceLength)
    ChimeraAverageScore=sum(ChimeraScore[(Boundary1-1):(SpliceLength+Boundary1)])/(SpliceLength)
    ScorePercentDifference=(OriginalAverageScore-ChimeraAverageScore)/OriginalAverageScore
    return ScorePercentDifference*100
def OverallConfidenceComparison(Protein):
    SARS2Score = list(map(float, open('SARS2.plddt', 'r').readlines()))
    ChimeraScore = list(map(float, open('SARS2w' + Protein + 'RBD.plddt', 'r').readlines()))
    SARS2AverageScore = sum(SARS2Score) / len(SARS2Score)
    ChimeraAverageScore = sum(ChimeraScore) / len(ChimeraScore)
    ScorePercentDifference = (SARS2AverageScore - ChimeraAverageScore) / SARS2AverageScore
    return ScorePercentDifference * 100

