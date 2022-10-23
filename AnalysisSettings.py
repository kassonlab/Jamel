import numpy as np
import Analysis
from RBDFinder import AlignmentFinder
import os
PresetList='Yes'
if PresetList=='Yes':
    BasenameList=[line.split()[-1] for line in open('List','r').readlines()]
    ProteinList=['3mer'+x.split()[-1] for x in BasenameList]
    AlignmentFileNames=[x+'onSARS2.aln' for x in BasenameList]
    PDB=[x+'.pdb' for x in ProteinList]
    #Plddtfiles=[x+'.plddt' for x in ProteinList]+ ['3merSARS2w'+x+'S1.plddt' for x in BasenameList]
    #PlddtResults=list(map(Analysis.AveragingMultimerPLDDT,Plddtfiles))


elif PresetList=='No':
    Plddtfiles=[x for x in os.listdir('/scratch/jws6pq/Notebook/Plddt/') if x[0]=='3']
    ProteinList=[x.replace('.plddt','') for x in Plddtfiles if x.find('3merSARS')==-1]
    Plddtfiles=list(map(Analysis.AveragingMultimerPLDDT,Plddtfiles))
    Plddtfiles=[x for x in Plddtfiles if x.find('3merSARS')==-1]
    BasenameList=[x.replace('3mer','') for x in ProteinList]
    AlignmentFileNames=[x.replace('3mer','')+'onSARS2.aln' for x in ProteinList if x.find('3merSARS')==-1]

os.chdir('/scratch/jws6pq/Notebook/Overall')
SequenceofInterest=['AYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLSEFRVYSSANNCTFEYVSQPFLKNLREFVFKNIDGYFKIYSKHTPPQGFSALEPLVDLPIGINITRFQTLLAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSYNYLYRNLKPFERDISTEIYNCYFPLQSYGFQPTVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICA' for x in ProteinList]
#S1= 1-540 AYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLSEFRVYSSANNCTFEYVSQPFLKNLREFVFKNIDGYFKIYSKHTPPQGFSALEPLVDLPIGINITRFQTLLAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSYNYLYRNLKPFERDISTEIYNCYFPLQSYGFQPTVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICA
#RBD= 224-425 TSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSYNYLYRNLKPFERDISTEIYNCYFPLQSYGFQPTVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVN
SARS2Splice1=[1 for x in ProteinList]
SARS2Splice2=[540 for x in ProteinList]
DomainSetting=['S1' for x in ProteinList]
ComparisonSetting=['3merSARS2' for x in ProteinList]

SpliceBoundaries=list(map(AlignmentFinder,AlignmentFileNames,SequenceofInterest))
ChimeraBoundary1=[x[0] for x in SpliceBoundaries]
ChimeraBoundary2=[x[1] for x in SpliceBoundaries]
OverallConfidence=list(map(Analysis.MultimerConfidenceComparison,ProteinList,ChimeraBoundary1,ChimeraBoundary2,SARS2Splice1,SARS2Splice2,DomainSetting,ComparisonSetting))
EmbossScore=list(map(Analysis.SequenceSimilarity,BasenameList,DomainSetting))
Overlap=Analysis.ContactOverlap,'SARS2wEverythingstable.aln'
FaultScan=list(map(Analysis.FaultScan,PDB))
DataChart=np.empty((len(ProteinList)+1,7),dtype=object)
DataChart[0,0]='Protein'
DataChart[0,1]='Sequence Similarity (%)'
DataChart[0,2]='Average Stability Difference'
DataChart[0,3]='Protein'
DataChart[0,4]='Overlap'
DataChart[0,5]='FullSimilarity'
DataChart[0,6]='PoorRender?'
DataChart[1:,0]=ProteinList
DataChart[1:,1]=EmbossScore
DataChart[1:,2]=OverallConfidence
DataChart[1:,3:6]=Overlap
DataChart[1:,6]=FaultScan
np.savetxt('/gpfs/gpfs0/scratch/jws6pq/CMfiles/'+DomainSetting[0]+'NewestChimeraAnalysis.tsv', DataChart, fmt="%s,%s,%s,%s,%s,%s,%s", delimiter="")
