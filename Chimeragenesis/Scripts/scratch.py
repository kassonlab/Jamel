import re

from AccessiontoAlignment import create_dictionary_from_alignment, create_seq_records, fasta_creation
from ChimeraGenerator import get_chimera_sequence
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
cuts=[152,350,416,647,821,996,1133,1274,1435] # RAndom CHimeras
# seq_dict=create_dictionary_from_alignment('..\Data\PDTvPDK_chimera.aln')
# dna_aln=create_dictionary_from_alignment('..\Data\PDKvPDT.aln')
# chimeras=dna_aln.copy()
# for x in cuts:
#     chimeras[f'PDTwPDK_{round(x/3)}']=chimeras['PDTH-2'].replace(chimeras['PDTH-2'][0:x],chimeras['PDK6'][0:x])
#     # print(chimeras['PDTH-2'].replace(chimeras['PDTH-2'][0:x],chimeras['PDK6'][0:x]))
# fasta_creation(f'PDTvPDK_dna.aln',sum((create_seq_records(label,seq) for label,seq in chimeras.items()),[]))
#
# dna_aln=create_dictionary_from_alignment('..\Data\PDTvPDK_dna.aln')
# for label,seq in dna_aln.items():
#
#     protein_seq=Seq(seq.replace('-','')).translate()
#     fasta_creation(f'PDTwPDK_{c}')
#     dna_aln[label]=protein_seq
#
# fasta_creation(f'../Data/PDTvPDK_chimera.aln', sum((create_seq_records(label, seq) for label,seq in dna_aln.items()), []))
# seq_dict=create_dictionary_from_alignment('..\Data\PDTvPDK_chimera.aln')
# labels=[]
# for label, seq in seq_dict.items():
#     cut=re.search(r'\d{2,9}',label)
#     if cut:
#         protein_cut=round(int(cut[0])/3)
#         label=re.sub(cut[0],str(protein_cut),label)
#     labels.append(f'/scratch/jws6pq/Notebook/ESM/RandomChimeras/{label}.fa')
#     fasta_creation(f'..\Data\\{label}.fa',create_seq_records(label,seq.replace('-',''),))
# for x in cuts:
#     cut = round(x / 3)
#     print(x,' ',f"{{'PDK6':(0,{cut}),'PDTH-2':({cut},-1)}}")
    # fasta_creation(f'..\Data\PDT0_{cut}.fa',
    #                create_seq_records(f'PDT0_{cut}', seq_dict['PDTH-2'].replace('-', '')[0:cut]))
# print(','.join(labels))
#
# alignment = PairwiseAligner()
# print(alignment.score(seq_dict['PDK6'], seq_dict['PDTH-2']) / alignment.score(seq_dict['PDK6'], seq_dict['PDK6']))

