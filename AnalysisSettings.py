import numpy as np
import Analysis
import RBDFinder
ProteinList=[line.split()[-1] for line in open('List','r').readlines()]
SARS2Splice1=['1' for x in ProteinList]
SARS2Splice2=['540' for x in ProteinList]
DomainSetting=['S1' for x in ProteinList]
AlignmentFileNames=[x+'onSARS2.aln' for x in ProteinList]
SequenceofInterest=['AYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLSEFRVYSSANNCTFEYVSQPFLKNLREFVFKNIDGYFKIYSKHTPPQGFSALEPLVDLPIGINITRFQTLLAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSYNYLYRNLKPFERDISTEIYNCYFPLQSYGFQPTVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICA' for x in ProteinList]
SpliceBoundaries=list(map(RBDFinder.AlignmentFinder,AlignmentFileNames,ProteinList,DomainSetting,SequenceofInterest))
ChimeraBoundary1=[x[0] for x in SpliceBoundaries]
ChimeraBoundary2=[x[1] for x in SpliceBoundaries]
RMSD=list(map(Analysis.PieceWiseRMSD,ProteinList,SARS2Splice1,SARS2Splice2,ChimeraBoundary1,ChimeraBoundary2))
DomainRMSD=[x[0] for x in RMSD]
OverallRMSD=[x[1] for x in RMSD]
PercentDifference=list(map(Analysis.SpliceConfidenceComparison,ProteinList,ChimeraBoundary1,ChimeraBoundary2,SARS2Splice1,DomainSetting))
OverallPercentDifference=list(map(Analysis.OverallConfidenceComparison,ProteinList,DomainSetting))
EmbossScore=list(map(Analysis.SequenceSimilarity,ProteinList,DomainSetting))
FoldXScore=list(map(Analysis.FoldXStability,ProteinList,DomainSetting))
DataChart=np.empty((len(ProteinList)+1,7),dtype=object)
DataChart[0,0]='Protein'
DataChart[0,1]=DomainSetting[0]+'RMSD'
DataChart[0,2]='OverallRMSD'
DataChart[0,3]='Emboss'
DataChart[0,4]=DomainSetting[0]+'Confidence'
DataChart[0,5]='OverallConfidence'
DataChart[0,6]='FoldX'
DataChart[1:,0]=ProteinList
DataChart[1:,1]=DomainRMSD
DataChart[1:,2]=OverallRMSD
DataChart[1:,3]=EmbossScore
DataChart[1:,4]=PercentDifference
DataChart[1:,5]=OverallPercentDifference
DataChart[1:,6]=FoldXScore
np.savetxt('/sfs/lustre/bahamut/scratch/jws6pq/CMfiles/'+DomainSetting[0]+'ChimeraAnalysis.tsv', DataChart, fmt="%s,%s,%s,%s,%s,%s,%s", delimiter="")
