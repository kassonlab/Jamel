import numpy as np
import Analysis
from AlignmentFinder import alignment_finder
import os
import concurrent.futures

basename_list=[line.split()[-1] for line in open('List', 'r').readlines()]
protein_list=[f'3mer{x.split()[-1]}' for x in basename_list]
alignment_files=[f'{protein}onSARS2.aln' for protein in basename_list]
pdb_files= [f'{protein}.pdb' for protein in protein_list] + [f'3merSARS2w{protein}S1.pdb' for protein in basename_list]
native_plddts = [f'Avg{protein}.plddt' for protein in protein_list]
chimera_plddt=[f'Avg3merSARS2w{protein}S1.plddt' for protein in basename_list]
plddt_files= [f'{protein}.plddt' for protein in protein_list] + [f'3merSARS2w{protein}S1.plddt' for protein in basename_list]

os.chdir('/scratch/jws6pq/Notebook/Overall')
SequenceofInterest=['AYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCND' \
                    'PFLSEFRVYSSANNCTFEYVSQPFLKNLREFVFKNIDGYFKIYSKHTPPQGFSALEPLVDLPIGINITRFQTLLAAYYVGYLQPRTFLLKYNENGTI' \
                    'TDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVS' \
                    'PTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSYNYLYRNLKPFERDISTEIYNCYFPLQSYGFQPTVGYQ' \
                    'PYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNT' \
                    'SNQVAVLYQDVNCTEVTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICA' for x in protein_list]
ComparisonSetting=['3merSARS2' for x in protein_list]
# PlddtResults=list(map(Analysis.AveragingMultimerPLDDT,Plddtfiles))

with concurrent.futures.ProcessPoolExecutor() as executor:
    SpliceBoundaries = list(executor.map(alignment_finder, alignment_files, SequenceofInterest))
    SpliceLength=[x[1]-x[0] for x in SpliceBoundaries]
    Similarity=list(executor.map(Analysis.get_sequence_similarity, [f'{protein}.emboss' for protein in basename_list]))
    AverageDifference=[] ; sars_difference=[] ; native_difference=[]
    for i in range(len(native_plddts)):
        Section1=Analysis.relative_stability('Avg3merSARS2.plddt', chimera_plddt[i], (0, 1), (0, 1))
        Section2=Analysis.relative_stability(native_plddts[i], chimera_plddt[i], (1, 1 + SpliceLength[i]),
                                             (SpliceBoundaries[i][0], SpliceBoundaries[i][1]))
        native_difference.append(Section2[0]/Section2[1])
        Section3=Analysis.relative_stability('Avg3merSARS2.plddt', chimera_plddt[i], (1 + SpliceLength[i], None), (540, None))
        sars_difference.append(Section3[0] / Section3[1])
        AverageRelativeDifference=(Section1[0]+Section2[0]+Section3[0])/(Section1[1]+Section2[1]+Section3[1])
        AverageDifference.append(AverageRelativeDifference)
    OverallDiff=list(executor.map(Analysis.overall_confidence, native_plddts))
    OverallChiDiff = list(executor.map(Analysis.overall_confidence, chimera_plddt))
DataChart=np.empty((len(protein_list) + 1, 7), dtype=object)
DataChart[0,0],DataChart[1:,0]='Protein', protein_list
DataChart[0,1],DataChart[1:,1]='S1 Sequence Similarity (%)',Similarity
DataChart[0,2],DataChart[1:,2]='Overall native plddt',OverallDiff
DataChart[0,3],DataChart[1:,3]='Overall chimera plddt',OverallChiDiff
DataChart[0,4],DataChart[1:,4]='Average Stability Difference',AverageDifference
DataChart[0,5],DataChart[1:,5]='SARS Difference',sars_difference
DataChart[0,6],DataChart[1:,6]='Native Difference',native_difference
np.savetxt(f'/gpfs/gpfs0/scratch/jws6pq/CMfiles/S1_0103ChimeraAnalysis.tsv', DataChart, fmt="%s,%s,%s,%s,%s,%s,%s", delimiter="")
