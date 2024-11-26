import ast
import os
import re
from enum import Enum
from json import load
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import torch
import umap
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from AccessiontoAlignment import create_dictionary_from_alignment, dictionary_to_fasta, \
    translate_dna_to_protein
from ChimeraGenerator import sequence_splice, chimera_sequence_creation, \
    create_chimera_combinations


class EmbedDims(Enum):
    Dim1 = 'mean'
    Dim2 = 'per_tok'

    def __str__(self):
        return self.name


def label_to_file(pt_directory, label):
    return os.path.join(pt_directory, f'{label}.pt')


def tensor_distance(tensor1, tensor2):
    dim = len(tuple(dim for dim in tensor1.shape if dim != 1))
    dist_func = torch.linalg.matrix_norm if dim == 2 else torch.linalg.vector_norm
    return abs(dist_func(tensor1[0] - tensor2[0])).item()


class SequenceDataframe(pd.DataFrame):
    def __init__(self):
        super().__init__(columns=['aln_sequence', 'sequence', 'description'])

    def dataframe_to_fa(self, new_fasta_file):
        seq_records = [
            SeqRecord(Seq(seq_info['aln_sequence']), description=str(seq_info['description']), id=str(seq_id))
            for seq_id, seq_info in self.iterrows()]
        with open(new_fasta_file, "w") as output_handle:
            SeqIO.write(seq_records, output_handle, "fasta")

    def create_dataframe_from_alignment(self, alignment_file):
        """Takes a fasta style alignment and makes a dictionary where the key is whatever signifier follows '>'
        and the value is the sequence with no spaces"""
        with open(alignment_file) as handle:
            for seq in SeqIO.parse(handle, "fasta"):
                self.loc[seq.id] = {'sequence': str(seq.seq).replace('-', ''), 'aln_sequence': str(seq.seq),
                                    'description': seq.description.replace(f'{seq.id} ', '')}
        return self

    def create_column(self, column_name):
        self[column_name] = np.nan

    def add_value(self, label, column, value):
        self.loc[label, column] = value

    def get_parents(self, chimera_label):
        return ast.literal_eval(self.loc[chimera_label, 'description'])

    def get_sequence(self, chimera_label):
        return self.loc[chimera_label, 'sequence']

    def get_aln(self, chimera_label):
        return self.loc[chimera_label, 'aln_sequence']


# with open(arg_json, 'rb') as jfile:
#     arg_dict = load(jfile)
# for key, value in arg_dict.items():
#     setattr(self, key, value)
#     def get_esm_embeddings(self):
#         if not os.path.exists(self.esm_output_directory):
#             os.makedirs(self.esm_output_directory)
#         os.system(
#             rf'python3 {self.extract_py_file} {self.esm_model} {self.chimera_seq_fasta}  {self.esm_output_directory} --include {self.embeddings_dimension.value}')
#         self.create_embedding_container()
# embeddings_dimension = 1
#     '''The ESM embeddings can be made per residue (2-dimensional) or just per sequence (1-dimensional)'''
# parent_aln_file: str
#     '''An fasta alignment file with your two parent proteins'''
#     splice_size: int
#     '''The size of the splice site that will move down the alignment at some step size x and recombine the two parent sequences'''
#     chimera_seq_fasta: str
#     '''A single fasta file that contains all the chimeric sequences created by some splice and step size. (Including the parents)'''
#     step_size: int = 1
#     '''How much the splice site steps over after a recombination'''
#  extract_py_file: str
#     '''Path to the python script file from esm scripts called extract.py'''
# esm_output_directory: str
#     '''Directory in create all esm outputs (.pt files) and the chimera_seq_output. If the directory doesn't exist it'll be created for you.'''
# esm_model: str = 'esm2_t30_150M_UR50D'
#     '''There are multiple different esm models with increasing amount of network layers and parameters. This
#     prediction models uses 150M parameter model for greatest performance (esm2_t48_15B_UR50D,esm2_t36_3B_UR50D,
#     esm2_t33_650M_UR50D,esm2_t30_150M_UR50D,esm2_t12_35M_UR50D,esm2_t6_8M_UR50D)'''
class EmbeddingAnalysis:
    distance_score_csv: str
    '''Csv file that contains the embedding distances'''
    EMBEDDINGS_DICT: dict[str:torch.Tensor]

    def __init__(self, embed_dict_pkl: str, aln_file):
        self.EMBEDDINGS_DICT: dict = torch.load(embed_dict_pkl)
        self.aln_df = SequenceDataframe()
        self.aln_df.create_dataframe_from_alignment(aln_file)

    def score_embedding_distance(self, chi_label, parent_labels):
        distance = 0
        for parent in parent_labels:
            distance += tensor_distance(self.EMBEDDINGS_DICT[chi_label], self.EMBEDDINGS_DICT[parent])
        if self.aln_df.get('distance') is not None and not self.aln_df.empty:
            self.aln_df.loc[chi_label, 'distance'] = distance
        return distance

    def per_residue_distance_for_val_set(self, chi_label, parent_labels):
        if self.EMBEDDINGS_DICT[chi_label].shape[0] != len(self.aln_df.loc[chi_label, 'sequence']):
            return
        distance = 0

        for parent, aln_positions in residue_inheritance.items():
            for pos in aln_positions:
                distance += tensor_distance(self.EMBEDDINGS_DICT[chi_label][pos, :],
                                            self.EMBEDDINGS_DICT[parent][pos, :])
        if self.aln_df.get('distance') is not None:
            self.aln_df.loc[chi_label, 'distance'] = distance
        return distance

    def score_all_per_res(self):
        self.aln_df['distance'] = np.nan
        for label, chi_tensor in self.EMBEDDINGS_DICT.items():
            self.per_residue_distance_for_val_set(label, self.aln_df.get_parents(label))

    def score_all_embeddings(self):
        self.aln_df['distance'] = np.nan
        for label, chi_tensor in self.EMBEDDINGS_DICT.items():
            self.score_embedding_distance(label, self.aln_df.get_parents(label))


def pt_to_tensor(pt_file, dim1_or_dim2: EmbedDims = EmbedDims.Dim1, model='esm2_t30_150M_UR50D') -> torch.Tensor:
    dim1_or_dim2 = dim1_or_dim2.value
    repr_layer = int(re.search(r't(\d*)_', model).group(1))
    if dim1_or_dim2 == 'mean':
        tensor = torch.load(pt_file)['mean_representations'][repr_layer]
    else:
        tensor = torch.load(pt_file)['representations'][repr_layer]
        tensor = add_row_of_zeroes(tensor)
    return tensor


def add_row_of_zeroes(tensor: torch.Tensor):
    return torch.cat((torch.zeros(1, tensor.size(1)), tensor), dim=0)


def pickle_embed_dict_schema(directory, new_pkl_file):
    embed_dict = {}
    for file in Path.iterdir(directory):
        embed_dict[file.stem] = torch.load(file)
    torch.save(embed_dict, new_pkl_file)


def embedding_umap(embed_pkl_file: str):
    reducer = umap.UMAP(n_components=2)
    embed_dict: dict[str, torch.Tensor] = torch.load(embed_pkl_file)
    # rows in umap matrix are samples/proteins
    embedding_matrix = np.vstack(tuple(embed_dict.values()))
    umap_vectors = reducer.fit_transform(embedding_matrix)
    plt.scatter(umap_vectors[:, 0], umap_vectors[:, 1], cmap='viridis', marker='+')
    return plt.gcf()


if __name__ == '__main__':
    schema_data = pd.read_csv(r"C:\Users\jamel\Downloads\schema_data.csv", index_col='chimera_block_ID')
    test = EmbeddingAnalysis(r'large_ankh.pkl', r'labeled_schema_aln')
    test.score_all_per_res()
    # pd.concat((schema_data, test.aln_df), axis=1).to_csv('large_ankh_res.csv')
    # ankh = torch.load('large_ankh.pkl')
    #
    # prost={label:tensor.transpose(0,1) for label,tensor in prost.items()}
    # print(ankh['c0000000000'].shape)
    # prost=torch.load('prost.pkl')
    #
    # prost={label:tensor.transpose(0,1) for label,tensor in prost.items()}
    # print(prost['c0000000000'].shape)
    # torch.save(prost,'prost.pkl')
    # plt.title('Riff CsChrim UMAP')
    # plt.xlabel('Component 1')
    # plt.ylabel('Component 2')
    # plt.colorbar(label='Embedding Distance')
    # plt.show()
    # parser = argparse.ArgumentParser(
    #     description='Creating a fasta file of all potential chimeric proteins between two parents based on a sliding splice site')
    # parser.add_argument('-in', '--inputjson', type=str, required=True,
    #                     help='alignment file in fasta style between 2 proteins')
    # args = parser.parse_args()
    # if args.inputjson:
    #     esm_args = EmbeddingAnalysis(args.inputjson)
    #     parent_seq_dict = create_dictionary_from_alignment(esm_args.parent_aln_file)
    #     chi_combo_dict = create_chimera_combinations(parent_seq_dict, esm_args.splice_size,
    #                                                  scanner_rate=esm_args.step_size,
    #                                                  new_fasta_file=esm_args.chimera_seq_fasta)
    #     if os.path.exists(esm_args.chimera_seq_fasta):
    #         # esm_args.get_esm_embeddings()
    #         esm_args.create_embedding_container()
    #         embeddings = esm_args.score_all_embeddings()
    #         scores = score_all_embeddings(esm_args.esm_output_directory, esm_args.chimera_seq_output)
    #         scores.to_csv(esm_args.distance_score_csv)
