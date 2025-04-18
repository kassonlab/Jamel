from collections import Counter
from pathlib import Path

import pandas as pd
import torch
import umap
from matplotlib import pyplot as plt, colormaps
from AccessiontoAlignment import create_dictionary_from_alignment, create_seq_records, fasta_creation, alignment_finder
from Bio.Seq import Seq

from ChimeraClasses import HomomerChimeraArgs
from ESM import SequenceDataframe
import numpy as np

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

# Looking at the umap of all 11 proteins for each residue position
# for x in range(45,schema.EMBEDDINGS_DICT['PDTwPDK_647'].shape[0]):
#     labels=[label for label in schema.EMBEDDINGS_DICT.keys()]
#     res_stack=[res_tensor[x,:].shape for res_tensor in schema.EMBEDDINGS_DICT.values()]
#     embedding_matrix = np.vstack(res_stack)
#     umap_3d = umap.UMAP(n_components=2)
#     X_umap = umap_3d.fit_transform(embedding_matrix)
#
#     # Create 3D scatter plot
#     plt.scatter(X_umap[:, 0], X_umap[:, 1], marker='+')
#     for i in range(len(X_umap[:, 0])):
#         plt.text(X_umap[i, 0], X_umap[i, 1], s=f"{labels[i]}", fontsize=10, ha="right", va="bottom")

    # # Labels and title
    #
    # ax.set_xlabel("UMAP 1")
    # ax.set_ylabel("UMAP 2")
    # ax.set_zlabel("UMAP 3")

    # plt.show()
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

df = pd.read_excel(r"C:\Users\jamel\PycharmProjects\Jamel\Chimeragenesis\Data\eid_sars_v5.xlsx", sheet_name='sars_v5')
print(df.columns)
df.plot.scatter(x='base_len', y='dot_product_rank',s=1)
df.plot.scatter(x='chimera_len', y='dot_product',s=1)
plt.show()