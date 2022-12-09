import numpy as np
import Analysis

Proteinlist=[line.split()[-1] for line in open('/gpfs/gpfs0/scratch/jws6pq/BridCMfiles/Flu_CovidList','r').readlines()]
PDB=[x+'.pdb' for x in Proteinlist]
Plddtfiles=[x+'.plddt' for x in Proteinlist]

PlddtResults=list(map(Analysis.AveragingMultimerPLDDT,Plddtfiles))

SpliceLength1=200
Boundary1=[x for x in range(0,489-SpliceLength1,10)]
Protein1=['HA.fasta' for x in Boundary1]
Boundary2=[x+SpliceLength1 for x in Boundary1]

SpliceLength2=200
Boundary3=[x for x in range(540,973-SpliceLength2,20)]
Protein2=['SARS2.fasta' for x in Boundary3]
Boundary4=[x+SpliceLength1 for x in Boundary3]

AverageRelativeDifference=[]
HADifference=[]
SARSDifference=[]
k=0
for i in range(len(Boundary1)):
    for j in range(len(Boundary3)):
        print(PlddtResults[k], Boundary1[i], Boundary2[i], Boundary3[j], Boundary4[j],k)
        Protein1AverageDifference = Analysis.MultimerConfidenceComparison('AvgHA.plddt', PlddtResults[k], [0, 0 + SpliceLength1],[Boundary1[i], Boundary2[i]])
        Protein2AverageDifference = Analysis.MultimerConfidenceComparison('Avg3merSARS2.plddt', PlddtResults[k],[0 + SpliceLength1, None],[Boundary3[j], Boundary4[j]])
        AverageRelativeDifference.append((Protein1AverageDifference[0] + Protein2AverageDifference[0]) / ( Protein1AverageDifference[1] + Protein2AverageDifference[1]))
        HADifference.append(Protein1AverageDifference[0]/Protein1AverageDifference[1])
        SARSDifference.append(Protein2AverageDifference[0]/Protein2AverageDifference[1])
        k+=1
OverallDiff=list(map(Analysis.OverallConfidence,Plddtfiles))

DataChart=np.empty((len(Proteinlist)+1,5),dtype=object)
DataChart[0,0],DataChart[1:,0]='Protein',Proteinlist
DataChart[0,1],DataChart[1:,1]='Average Stability Difference',AverageRelativeDifference
DataChart[0,2],DataChart[1:,2]='HA Difference',HADifference
DataChart[0,3],DataChart[1:,3]='SARS Difference',SARSDifference
DataChart[0,4],DataChart[1:,4]='Overall chimera plddt',OverallDiff


np.savetxt('/gpfs/gpfs0/scratch/jws6pq/CMfiles/120522HACovidAnalysis.tsv', DataChart, fmt="%s,%s,%s,%s,%s", delimiter="")
