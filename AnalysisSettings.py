import numpy as np
import Analysis
import RBDFinder
ProteinList=[line.split()[-1] for line in open('List','r').readlines()]
SARS2Splice1=[224 for x in ProteinList]
SARS2Splice2=[425 for x in ProteinList]
DomainSetting=['RBD' for x in ProteinList]
AlignmentFileNames=[x+'onSARS2.aln' for x in ProteinList]
SequenceofInterest=['TSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSYNYLYRNLKPFERDISTEIYNCYFPLQSYGFQPTVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVN' for x in ProteinList]
#S1= 1-540 AYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLSEFRVYSSANNCTFEYVSQPFLKNLREFVFKNIDGYFKIYSKHTPPQGFSALEPLVDLPIGINITRFQTLLAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSYNYLYRNLKPFERDISTEIYNCYFPLQSYGFQPTVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICA
#RBD= 224-425 TSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSYNYLYRNLKPFERDISTEIYNCYFPLQSYGFQPTVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVN

SpliceBoundaries=list(map(RBDFinder.AlignmentFinder,AlignmentFileNames,ProteinList,DomainSetting,SequenceofInterest))
ChimeraBoundary1=[x[0] for x in SpliceBoundaries]
ChimeraBoundary2=[x[1] for x in SpliceBoundaries]
RMSD=list(map(Analysis.PieceWiseRMSD,ProteinList,SARS2Splice1,SARS2Splice2,ChimeraBoundary1,ChimeraBoundary2))
DomainRMSD=[x[0] for x in RMSD]
OverallRMSD=[x[1] for x in RMSD]
PercentDifference=list(map(Analysis.SpliceConfidenceComparison,ProteinList,ChimeraBoundary1,ChimeraBoundary2,SARS2Splice1,SARS2Splice2,DomainSetting))
OverallConfidence=list(map(Analysis.OverallConfidenceComparison,ProteinList,DomainSetting))
OverallDifference=[x[0] for x in OverallConfidence]
OverallScore=[x[1] for x in OverallConfidence]
EmbossScore=list(map(Analysis.SequenceSimilarity,ProteinList,DomainSetting))
FoldXScore=list(map(Analysis.FoldXStability,ProteinList,DomainSetting))
DataChart=np.empty((len(ProteinList)+1,8),dtype=object)
DataChart[0,0]='Protein'
DataChart[0,1]=DomainSetting[0]+'RMSD'
DataChart[0,2]='Overall RMSD'
DataChart[0,3]='Sequence Similarity (%)'
DataChart[0,4]=DomainSetting[0]+'Stability Difference (%)'
DataChart[0,5]='Overall Stability Difference'
DataChart[0,6]='Overall Stability'
DataChart[0,7]='FoldX'
DataChart[1:,0]=ProteinList
DataChart[1:,1]=DomainRMSD
DataChart[1:,2]=OverallRMSD
DataChart[1:,3]=EmbossScore
DataChart[1:,4]=PercentDifference
DataChart[1:,5]=OverallDifference
DataChart[1:,6]=OverallScore
DataChart[1:,7]=FoldXScore

np.savetxt('/sfs/lustre/bahamut/scratch/jws6pq/CMfiles/'+DomainSetting[0]+'NewChimeraAnalysis.tsv', DataChart, fmt="%s,%s,%s,%s,%s,%s,%s,%s", delimiter="")
