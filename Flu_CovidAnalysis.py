import numpy as np
import Analysis
from RBDFinder import AlignmentFinder
import os
import concurrent.futures
PresetList='Yes'
Proteinlist=[line.split()[-1] for line in open('Flu_CovidList','r').readlines()]
PDB=[x+'.pdb' for x in Proteinlist]
Plddtfiles=[x+'.plddt' for x in Proteinlist]

os.chdir('/scratch/jws6pq/Notebook/Overall')
PlddtResults=list(map(Analysis.AveragingMultimerPLDDT,Plddtfiles))

SpliceLength1=200
Boundary1=[x for x in range(0,489-SpliceLength1,10)]
Protein1=['HA.fasta' for x in Boundary1]
Boundary2=[x+SpliceLength1 for x in Boundary1]

SpliceLength2=200
Boundary3=[x for x in range(540,973-SpliceLength2,20)]
Protein2=['SARS2.fasta' for x in Boundary3]
Boundary4=[x+SpliceLength1 for x in Boundary3]

with concurrent.futures.ProcessPoolExecutor() as executor:
    AverageDifference=list(executor.map(Analysis.MultimerConfidenceComparison,Proteinlist,ChimeraBoundary1,ChimeraBoundary2,SARS2Splice1,SARS2Splice2,DomainSetting,ComparisonSetting))
    OverallDiff=list(executor.map(Analysis.OverallConfidence,Plddtfiles))
    # Overlap=list(executor.map(Analysis.ContactOverlap,OverlapAlignment,BasenameList))
OverlapScore=[x[1] for x in Overlap]
OverallSimilarity=[x[2] for x in Overlap]
DataChart=np.empty((len(ProteinList)+1,6),dtype=object)
DataChart[0,0],DataChart[1:,0]='Protein',ProteinList
DataChart[0,1],DataChart[1:,1]='Average Stability Difference',AverageDifference
DataChart[0,2],DataChart[1:,2]='Overall native plddt',OverallDiff[0:len(ProteinList)]
DataChart[0,3],DataChart[1:,3]='Overall chimera plddt',OverallDiff[len(ProteinList):]
DataChart[0,4],DataChart[1:,4]='Sequence Similarity (%)',OverallSimilarity
DataChart[0,5],DataChart[1:,5]='Contact Overlap',OverlapScore

np.savetxt('/gpfs/gpfs0/scratch/jws6pq/CMfiles/'+DomainSetting[0]+'1031ChimeraAnalysis.tsv', DataChart, fmt="%s,%s,%s,%s,%s,%s", delimiter="")
