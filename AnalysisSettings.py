import numpy as np
import Analysis
import RBDFinder
import os

PresetList='No'
if PresetList=='Yes':
    ProteinList=[line.split()[-1] for line in open('List1','r').readlines()]

elif PresetList=='No':
    Plddtfiles=os.listdir('/sfs/lustre/bahamut/scratch/jws6pq/Notebook/Plddt/')
    Plddtfiles=[x for x in Plddtfiles if x[0]=='3']
    ProteinList=[x.replace('.plddt','') for x in Plddtfiles if x.find('3merSARS')==-1]
    Plddtfiles=list(map(Analysis.AveragingMultimerPLDDT,Plddtfiles))
    Plddtfiles=[x for x in Plddtfiles if x.find('3merSARS')==-1]

SARS2Splice1=[1 for x in ProteinList]
SARS2Splice2=[540 for x in ProteinList]
DomainSetting=['S1' for x in ProteinList]
ComparisonSetting=['3merSARS2' for x in ProteinList]

os.chdir('/sfs/lustre/bahamut/scratch/jws6pq/Notebook/Overall')
AlignmentFileNames=[x.replace('3mer','')+'onSARS2.aln' for x in ProteinList if x.find('3merSARS')==-1]
SequenceofInterest=['AYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLSEFRVYSSANNCTFEYVSQPFLKNLREFVFKNIDGYFKIYSKHTPPQGFSALEPLVDLPIGINITRFQTLLAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSYNYLYRNLKPFERDISTEIYNCYFPLQSYGFQPTVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICA' for x in ProteinList]
#S1= 1-540 AYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLSEFRVYSSANNCTFEYVSQPFLKNLREFVFKNIDGYFKIYSKHTPPQGFSALEPLVDLPIGINITRFQTLLAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSYNYLYRNLKPFERDISTEIYNCYFPLQSYGFQPTVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICA
#RBD= 224-425 TSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSYNYLYRNLKPFERDISTEIYNCYFPLQSYGFQPTVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVN

SpliceBoundaries=list(map(RBDFinder.AlignmentFinder,AlignmentFileNames,SequenceofInterest))
ChimeraBoundary1=[x[0] for x in SpliceBoundaries]
ChimeraBoundary2=[x[1] for x in SpliceBoundaries]
BasenameList=[x.replace('3mer','') for x in ProteinList]
OverallConfidence=list(map(Analysis.MultimerConfidenceComparison,ProteinList,ChimeraBoundary1,ChimeraBoundary2,SARS2Splice1,SARS2Splice2,DomainSetting,ComparisonSetting))
EmbossScore=list(map(Analysis.SequenceSimilarity,BasenameList,DomainSetting))
FoldXScore=list(map(Analysis.FoldXStability,BasenameList,DomainSetting))
DataChart=np.empty((len(ProteinList)+1,4),dtype=object)
DataChart[0,0]='Protein'

DataChart[0,1]='Sequence Similarity (%)'
DataChart[0,2]='Average Stability Difference'
DataChart[0,3]='FoldX'
DataChart[1:,0]=ProteinList

DataChart[1:,1]=EmbossScore
DataChart[1:,2]=OverallConfidence
DataChart[1:,3]=FoldXScore
np.savetxt('/sfs/lustre/bahamut/scratch/jws6pq/CMfiles/'+DomainSetting[0]+'NewestChimeraAnalysis.tsv', DataChart, fmt="%s,%s,%s,%s", delimiter="")
