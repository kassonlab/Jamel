def ConfidenceComparison(Protein,Boundary1,Boundary2,Boundary3,Boundary4):
    import math
    import numpy as np
    ProteinScore=list(map(float,open(Protein+'.plddt', 'r').readlines()))
    SARS2Score=list(map(float,(open('SARS2.plddt', 'r').readlines())))
    ChimeraScore=list(map(float,open('SARS2w'+Protein+'RBD.plddt', 'r').readlines()))
    SpliceLength=Boundary4-Boundary3
    OriginalAverageScore=sum(ProteinScore[(Boundary3-1):(Boundary4)])/(SpliceLength)
    ChimeraAverageScore=sum(ChimeraScore[(Boundary1-1):(SpliceLength+Boundary1)])/(SpliceLength)
    ScorePercentDifference=(OriginalAverageScore-ChimeraAverageScore)/OriginalAverageScore
    return ScorePercentDifference*100
import RBDFinder
SpliceBoundary=RBDFinder.AlignmentFinder('RhinoAlphaonSARS2.aln','RhinoAlpha')
ConfidenceComparison('RhinoAlpha',234,432,SpliceBoundary[0],SpliceBoundary[1])
