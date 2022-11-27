import numpy as np
import Analysis
from RBDFinder import AlignmentFinder
import os
import concurrent.futures
PresetList='Yes'
if PresetList=='Yes':
    BasenameList=[line.split()[-1] for line in open('List','r').readlines()]
    ProteinList=['3mer'+x.split()[-1] for x in BasenameList]
    AlignmentFileNames=[x+'onSARS2.aln' for x in BasenameList]
    PDB=[x+'.pdb' for x in ProteinList]+['3merSARS2w'+x+'S1.pdb' for x in BasenameList]
    Plddtfiles = [('Avg'+x + '.plddt','Avg3merSARS2w' + x + 'S1.plddt') for x in ProteinList]
    # Plddtfiles=[x+'.plddt' for x in ProteinList]+['3merSARS2w'+x+'S1.plddt' for x in BasenameList]
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
#S1= 0-539 AYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLSEFRVYSSANNCTFEYVSQPFLKNLREFVFKNIDGYFKIYSKHTPPQGFSALEPLVDLPIGINITRFQTLLAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSYNYLYRNLKPFERDISTEIYNCYFPLQSYGFQPTVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICA
#RBD= 223-424 TSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSYNYLYRNLKPFERDISTEIYNCYFPLQSYGFQPTVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVN
SARS2Splice=[(1,540) for x in ProteinList]
DomainSetting=['S1' for x in ProteinList]
ComparisonSetting=['3merSARS2' for x in ProteinList]
OverlapAlignment=['SARS2wEverythingstable.aln' for x in BasenameList]


with concurrent.futures.ProcessPoolExecutor() as executor:
    SpliceBoundaries = list(executor.map(AlignmentFinder, AlignmentFileNames, SequenceofInterest))

    AverageDifference=list(executor.map(Analysis.MultimerConfidenceComparison,ProteinList,ChimeraBoundary1,ChimeraBoundary2,SARS2Splice1,SARS2Splice2,DomainSetting,ComparisonSetting))
    OverallDiff=list(executor.map(Analysis.OverallConfidence,Plddtfiles))
    Overlap=list(executor.map(Analysis.ContactOverlap,OverlapAlignment,BasenameList))
    FaultScan=list(executor.map(Analysis.FaultScan,PDB))
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
