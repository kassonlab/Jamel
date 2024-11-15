import ast
import os
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from matplotlib import pyplot as plt
import ESM
from AccessiontoAlignment import create_dictionary_from_alignment, dictionary_to_fasta
from Chimeragenesis.ESM import label_to_file, pt_to_tensor, EmbedDims

parent_ids = {'0': 'c0000000000', '1': 'c1111111111', '2': 'c2222222222'}
SCHEMA_ALN = r'C:\Users\jamel\PycharmProjects\Jamel\Chimeragenesis\esm\schema_msa.aln'
schema_data = pd.read_csv(r'C:\Users\jamel\PycharmProjects\Jamel\Chimeragenesis\schema_data.csv').set_index(
    'chimera_block_ID')
schema_data[['residue_dist_1d', 'residue_dist_2d']] = np.nan


class SchemaScorer:
    def __init__(self, chimera_label, aln_file, pt_direc, model):
        self.model = model
        self.distance = None
        self.chimera_label = chimera_label
        self.aln_file = aln_file
        self.seq_df = self.aln_to_pd()
        self.similarity_dict = {}
        self.residue_similarity_cipher()
        # self.embed_dict = self.aln_to_embedding(pt_direc)
        self.pt_direc = pt_direc

    def aln_to_pd(self):
        seq_dict = create_dictionary_from_alignment(self.aln_file)
        seq_dict = {label: list(seq) for label, seq in seq_dict.items()}
        return pd.DataFrame(seq_dict)

    def residue_similarity_cipher(self) -> dict[int, list]:
        """This takes an alignment dataframe from aln_to_pd and creates a dictionary that has alignment position as keys
        and a list of the parents that share amino acids at that position with the chimera/child"""
        parents = self.seq_df.columns.to_list()
        parents.remove(self.chimera_label)
        for row, position in self.seq_df.iterrows():
            self.similarity_dict[int(row)] = [parent for parent in parents if
                                              position[self.chimera_label] == position[parent]]
        return self.similarity_dict

    def score_distance(self,dims: ESM.EmbedDims):
        dist_func = torch.linalg.matrix_norm if dims == EmbedDims.Dim2 else torch.linalg.vector_norm
        seq_dict = create_dictionary_from_alignment(self.aln_file)
        embed_dict = {label: pt_to_tensor(label_to_file(self.pt_direc, label), dims, self.model) for label in
                      seq_dict.keys()}
        parent1, parent2 = (label for label in embed_dict if label != self.chimera_label)
        paren1_dist = abs(dist_func(embed_dict[self.chimera_label] - embed_dict[parent1]))
        paren2_dist = abs(dist_func(embed_dict[self.chimera_label] - embed_dict[parent2]))
        return paren1_dist + paren2_dist


def separate_aln():
    schema_data['parent_label'] = schema_data['parent_label'].apply(lambda label: ast.literal_eval(label))
    schema_seq_dict = create_dictionary_from_alignment(SCHEMA_ALN)
    for chi_id, data in schema_data.iterrows():
        if len(data['parent_label']) == 3:
            parent1, parent2 = (parent_ids[parent] for parent in data['parent_label'][:2])
            new_seq_dict = {chi_id: schema_seq_dict[chi_id], parent1: schema_seq_dict[parent1],
                            parent2: schema_seq_dict[parent2]}
            dictionary_to_fasta(new_seq_dict,
                                fr'C:\Users\jamel\PycharmProjects\Jamel\Chimeragenesis\schema_individs\{chi_id}.aln')


if __name__ == '__main__':
    # ESM.ESMArguments('schema.json').get_esm_embeddings()
    # separate_aln()
    model = 'esm2_t36_3B_UR50D'
    aln_direc = r'C:\Users\jamel\PycharmProjects\Jamel\Chimeragenesis\schema_individs'
    for aln in os.listdir(aln_direc):
        label = Path(aln).stem
        chimera = SchemaScorer(label, os.path.join(aln_direc, aln),
                               r'C:\Users\jamel\PycharmProjects\Jamel\esm\3B_model\1d_schema_data', model)
        twod_chimera = SchemaScorer(label, os.path.join(aln_direc, aln),
                                    r'C:\Users\jamel\PycharmProjects\Jamel\esm\3B_model\2d_schema_data', model)
        schema_data.loc[label, 'residue_dist_1d'] = chimera.score_distance(EmbedDims.Dim1).item()
        schema_data.loc[label, 'residue_dist_2d'] = twod_chimera.score_distance(EmbedDims.Dim2).item()
    schema_data.to_csv('3B.csv')
