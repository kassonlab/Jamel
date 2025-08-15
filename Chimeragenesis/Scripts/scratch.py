import ast
import re
from collections import Counter
from datetime import date
from pathlib import Path
import pandas as pd
import torch, umap, Analysis
from matplotlib import pyplot as plt
from AccessiontoAlignment import create_dictionary_from_alignment, create_seq_records, fasta_creation, alignment_finder,calculate_sequence_identity,convert_fasta_aln_to_clustal
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
from ChimeraClasses import HomomerChimeraArgs
from Analysis import get_plddt_dict_from_pdb
from ESM import SequenceDataframe
import numpy as np
from pandasgui import show

cuts=[152,350,416,647,821,996,1133,1274,1435] # RAndom CHimeras
# seq_dict=create_dictionary_from_alignment('..\Data\PDTvPDK_chimera.aln')
# dna_aln=create_dictionary_from_alignment('..\Data\PDKvPDT_dna.aln')
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

# aln_file='..\Data\PDTvPDK_chimera.aln'
# seq_df = SequenceDataframe()
# seq_dict = create_dictionary_from_alignment(aln_file)
# for x in cuts:
#     cut = round(x / 3)
#     label=f'PDTwPDK_{x}'
#     inheritance=alignment_finder(seq_dict['PDTH-2'].replace('-','')[0:cut], 'PDK6', 'PDTH-2', aln_file)[1]
#     for parent,bounds in inheritance.items():
#         for parent_bound,chi_bound in bounds.items():
#             print(parent)
#             print(seq_dict[parent].replace('-','')[slice(*parent_bound)])
#             print(seq_dict[label].replace('-','')[slice(*chi_bound)])

#     seq_df.add_protein(label, seq_dict[label], inheritance)
# seq_df.dataframe_to_aln(Path(aln_file).with_suffix('.inh'))
    # for label, seq in seq_dict.items():
#     print(x,' ',f"{{'PDK6':(0,{cut}),'PDTH-2':({cut},-1)}}")
    # fasta_creation(f'..\Data\PDT0_{cut}.fa',
    #                create_seq_records(f'PDT0_{cut}', seq_dict['PDTH-2'].replace('-', '')[0:cut]))
# print(','.join(labels))
#
# randoms=SequenceDataframe(r'/Chimeragenesis/Data/Randoms/PDTvPDK_chimera.inh', embed_dict_pkl=r'C:\Users\jamel\PycharmProjects\Jamel\Chimeragenesis\Data\randoms_prost_X.pkl')


# Creating inheritance file for schmea
# schema_cuts=[0,35,49,94,117,146,189,214,232,266,None]
# schema_dict=create_dictionary_from_alignment(r'C:\Users\jamel\PycharmProjects\Jamel\Chimeragenesis\Data\notag_schema_parents.aln')
# notag_df=SequenceDataframe(r'C:\Users\jamel\PycharmProjects\Jamel\Chimeragenesis\Data\schema_no_tag.aln')
# parent_dict={'0':'c0000000000','1':'c1111111111','2':'c2222222222'}
# for label in notag_df.index:
#     label=str(label).replace('c','')
#     counts=Counter(label)
#     partner_symbol=min(counts,key=counts.get)
#     base_symbol=max(counts,key=counts.get)
#     block_index=label.find(partner_symbol)
#     seq_of_int=schema_dict[parent_dict[base_symbol]][slice(*schema_cuts[block_index:block_index+2])].replace('-','')
#     _,inheritance=alignment_finder(seq_of_int,parent_dict[partner_symbol],parent_dict[base_symbol],r'C:\Users\jamel\PycharmProjects\Jamel\Chimeragenesis\Data\notag_schema_parents.aln')
#     notag_df.add_value('c'+label,'description',inheritance)
# notag_df.dataframe_to_aln(r'C:\Users\jamel\PycharmProjects\Jamel\Chimeragenesis\Data\notag_schema_parents.inh')

# aligner = PairwiseAligner()
# aligner.mode = 'global'
# fig, ax = plt.subplots()
# df = pd.read_excel(r"C:\Users\jamel\PycharmProjects\Jamel\Chimeragenesis\Data\eid_sars_v5.xlsx", sheet_name='sars',index_col='label')
# print(df.columns)
# # get sequence Identity
# df['normalized_dp']=df['dot_product']/df['chimera_len']
# for label,data in df.iterrows():
#     if label in ['6VSB', 'EidolonBat']: continue
#     identity=[]
#     for parent in [df.loc['6VSB','sequence'],df.loc['EidolonBat','sequence']]:
#         identity.append(calculate_sequence_identity(*aligner.align(parent, data['sequence'])[0]))
#     df.loc[label, 'identity']=max(identity)
#
# df.plot.scatter(x='identity', y='normalized_dp',s=1,ax=ax)
# # df.to_csv('../Data/sars_w_identity.csv')
# plt.show()
# for label,data in df.iterrows():
#     if label in ['6VSB', 'EidolonBat']: continue
#     identity=[]
#     for parent in [df.loc['6VSB','sequence'],df.loc['EidolonBat','sequence']]:
#         identity.append(calculate_sequence_identity(*aligner.align(parent, data['sequence'])[0]))
#     df.loc[label, 'identity']=max(identity)



# df =pd.read_csv(r"..\Data\tem_cmy_aligned_af.csv",index_col='label').sort_values(by='dot_product').iloc[::42]
# print(*df['dot_product'],sep='\n')
# print(*df.index,sep='\n')
# final_7=["CMY1",
# "TEM1",
# "CMY1_TEM1_142_299",
# "CMY1_TEM1_163_326",
# "CMY1_TEM1_12_160",
# "CMY1_TEM1_26_235",
# "CMY1_TEM1_138_307"]
# tpr="/sfs/weka/scratch/jws6pq/Notebook/ESM/B-lac/aligned_cmy_tem/CMY1_chis/Gromacs/{0}/{0}_prod.tpr"
# xtc="/sfs/weka/scratch/jws6pq/Notebook/ESM/B-lac/aligned_cmy_tem/CMY1_chis/Gromacs/{0}/{0}_center.xtc"
# xvg="/sfs/weka/scratch/jws6pq/Notebook/ESM/B-lac/aligned_cmy_tem/CMY1_chis/Gromacs/{0}/{0}.xvg"
# print({label:xvg.format(label) for label in final_7})
show(pd.read_csv(r"C:\Users\jamel\Downloads\01mar2025-31may2025_form13f\INFOTABLE.tsv",delimiter='\t'))
# print(pd.read_csv(r"C:\Users\jamel\Downloads\01mar2025-31may2025_form13f\INFOTABLE.tsv"))

