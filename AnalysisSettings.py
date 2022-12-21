import numpy as np
import Analysis
from RBDFinder import alignment_finder
import os
import concurrent.futures
PresetList='Yes'
if PresetList=='Yes':
    BasenameList=[line.split()[-1] for line in open('List','r').readlines()]
    ProteinList=['3mer'+x.split()[-1] for x in BasenameList]
    AlignmentFileNames=[x+'onSARS2.aln' for x in BasenameList]
    PDB=[x+'.pdb' for x in ProteinList]+['3merSARS2w'+x+'S1.pdb' for x in BasenameList]
    NativePlddtfiles = ['Avg'+x + '.plddt' for x in ProteinList]
    ChimeraPlddtfiles=['Avg3merSARS2w' + x + 'S1.plddt' for x in BasenameList]
    Plddtfiles=[x+'.plddt' for x in ProteinList]+['3merSARS2w'+x+'S1.plddt' for x in BasenameList]



# elif PresetList=='No':
#     Plddtfiles=[x for x in os.listdir('/scratch/jws6pq/Notebook/Plddt/') if x[0]=='3']
#     ProteinList=[x.replace('.plddt','') for x in Plddtfiles if x.find('3merSARS')==-1]
#     Plddtfiles=list(map(Analysis.AveragingMultimerPLDDT,Plddtfiles))
#     Plddtfiles=[x for x in Plddtfiles if x.find('3merSARS')==-1]
#     BasenameList=[x.replace('3mer','') for x in ProteinList]
#     AlignmentFileNames=[x.replace('3mer','')+'onSARS2.aln' for x in ProteinList if x.find('3merSARS')==-1]
os.chdir('/scratch/jws6pq/Notebook/Overall')
SequenceofInterest=['AYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLSEFRVYSSANNCTFEYVSQPFLKNLREFVFKNIDGYFKIYSKHTPPQGFSALEPLVDLPIGINITRFQTLLAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSYNYLYRNLKPFERDISTEIYNCYFPLQSYGFQPTVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICA' for x in ProteinList]
Chimerasplice=[(0,539) for x in ProteinList]
DomainSetting=['S1' for x in ProteinList]
ComparisonSetting=['3merSARS2' for x in ProteinList]
# PlddtResults=list(map(Analysis.AveragingMultimerPLDDT,Plddtfiles))

with concurrent.futures.ProcessPoolExecutor() as executor:
    SpliceBoundaries = list(executor.map(alignment_finder, AlignmentFileNames, SequenceofInterest))
    SpliceLength=[x[1]-x[0] for x in SpliceBoundaries]
    Similarity=list(executor.map(Analysis.SequenceSimilarity, BasenameList, DomainSetting))
    AverageDifference=[]
    i=0
    for x in NativePlddtfiles:
        Section1=Analysis.MultimerConfidenceComparison('Avg3merSARS2.plddt',ChimeraPlddtfiles[i],(0,1),(0,1))
        Section2=Analysis.MultimerConfidenceComparison(NativePlddtfiles[i],ChimeraPlddtfiles[i],(1,1+SpliceLength[i]),(SpliceBoundaries[i][0],SpliceBoundaries[i][1]))
        Section3=Analysis.MultimerConfidenceComparison('Avg3merSARS2.plddt', ChimeraPlddtfiles[i], (1+SpliceLength[i],None), (540,None))
        AverageRelativeDifference=(Section1[0]+Section2[0]+Section3[0])/(Section1[1]+Section2[1]+Section3[1])
        AverageDifference.append(AverageRelativeDifference)
        i+=1
    OverallDiff=list(executor.map(Analysis.OverallConfidence,NativePlddtfiles))
    OverallChiDiff = list(executor.map(Analysis.OverallConfidence, ChimeraPlddtfiles))
DataChart=np.empty((len(ProteinList)+1,5),dtype=object)
DataChart[0,0],DataChart[1:,0]='Protein',ProteinList
DataChart[0,1],DataChart[1:,1]='Average Stability Difference',AverageDifference
DataChart[0,2],DataChart[1:,2]='Overall native plddt',OverallDiff
DataChart[0,3],DataChart[1:,3]='Overall chimera plddt',OverallChiDiff
DataChart[0,4],DataChart[1:,4]='S1 Sequence Similarity (%)',Similarity

np.savetxt('/gpfs/gpfs0/scratch/jws6pq/CMfiles/'+DomainSetting[0]+'_1101ChimeraAnalysis.tsv', DataChart, fmt="%s,%s,%s,%s,%s", delimiter="")
