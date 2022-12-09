import numpy as np
import Analysis
from RBDFinder import AlignmentFinder
import concurrent.futures

Listtoanalyze='List'
BasenameList=[line.split()[-1] for line in open(Listtoanalyze,'r').readlines()]
ProteinList=['3mer'+x.split()[-1] for x in BasenameList]
AlignmentFileNames=[x+'onSARS2.aln' for x in BasenameList]
PDB=[x+'.pdb' for x in ProteinList]+['3merSARS2w'+x+'RBD.pdb' for x in BasenameList]
NativePlddtfiles = ['Avg'+x + '.plddt' for x in ProteinList]
ChimeraPlddtfiles=['Avg3merSARS2w' + x + 'RBD.plddt' for x in BasenameList]
Plddtfiles=[x+'.plddt' for x in ProteinList]+['3merSARS2w'+x+'RBD.plddt' for x in BasenameList]
SequenceofInterest=['TSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSYNYLYRNLKPFERDISTEIYNCYFPLQSYGFQPTVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVN' for x in ProteinList]
Chimerasplice=[(223,424) for x in ProteinList]
DomainSetting=['RBD' for x in ProteinList]
ComparisonSetting=['3merSARS2' for x in ProteinList]
PlddtResults=list(map(Analysis.AveragingMultimerPLDDT,Plddtfiles))

with concurrent.futures.ProcessPoolExecutor() as executor:
    SpliceBoundaries = list(executor.map(AlignmentFinder, AlignmentFileNames, SequenceofInterest))
    SpliceLength=[x[1]-x[0] for x in SpliceBoundaries]
    Similarity=list(executor.map(Analysis.SequenceSimilarity, BasenameList, DomainSetting))
    AverageDifference=[]
    Fault=list(executor.map(Analysis.FaultScan,PDB))
    i=0
    for x in NativePlddtfiles:
        Section1=Analysis.MultimerConfidenceComparison('Avg3merSARS2.plddt',ChimeraPlddtfiles[i],(0,223),(0,223))
        Section2=Analysis.MultimerConfidenceComparison(NativePlddtfiles[i],ChimeraPlddtfiles[i],(223,223+SpliceLength[i]),(SpliceBoundaries[i][0],SpliceBoundaries[i][1]))
        Section3=Analysis.MultimerConfidenceComparison('Avg3merSARS2.plddt', ChimeraPlddtfiles[i], (223+SpliceLength[i],None), (424,None))
        AverageRelativeDifference=(Section1[0]+Section2[0]+Section3[0])/(Section1[1]+Section2[1]+Section3[1])
        AverageDifference.append(AverageRelativeDifference)
        i+=1
    OverallDiff=list(executor.map(Analysis.OverallConfidence,NativePlddtfiles))
    OverallChiDiff = list(executor.map(Analysis.OverallConfidence, ChimeraPlddtfiles))
DataChart=np.empty((len(ProteinList)+1,6),dtype=object)
DataChart[0,0],DataChart[1:,0]='Protein',ProteinList
DataChart[0,1],DataChart[1:,1]='Average Stability Difference',AverageDifference
DataChart[0,2],DataChart[1:,2]='Overall native plddt',OverallDiff
DataChart[0,3],DataChart[1:,3]='Overall chimera plddt',OverallChiDiff
DataChart[0,4],DataChart[1:,4]='RBD Sequence Similarity (%)',Similarity
DataChart[0,5],DataChart[1:,5]='FaultNative',Fault[0:len(ProteinList)]
DataChart[0,5],DataChart[1:,5]='FaultChi',Fault[len(ProteinList):]
np.savetxt('/gpfs/gpfs0/scratch/jws6pq/CMfiles/'+DomainSetting[0]+'_1201ChimeraAnalysis.tsv', DataChart, fmt="%s,%s,%s,%s,%s,%s,%s", delimiter="")