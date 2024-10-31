import os
import re
from enum import Enum
from json import load

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import torch
import umap
import argparse
from AccessiontoAlignment import create_dictionary_from_alignment, dictionary_to_fasta, \
    translate_dna_to_protein
from ChimeraGenerator import sequence_splice, chimera_sequence_creation, \
    create_chimera_combinations

# pt=torch.load(r'C:\Users\jamel\PycharmProjects\Jamel\Chimeragenesis\riffChr_data\Chrim_test.pt')
    # pt['parent1']=1
    # pt['parent2']=2
    # pt['splice_location']=3
    # print(pt)
    # torch.save(pt,r'C:\Users\jamel\PycharmProjects\Jamel\Chimeragenesis\riffChr_data\Chrim_test.pt')
class EmbedDims(Enum):
    Dim1 = 'mean'
    Dim2 = 'per_tok'

    def __str__(self):
        return self.name


def pt_to_tensor(pt_file, dim1_or_dim2: EmbedDims, repr_layer: int = 30) -> torch.Tensor:
    dim1_or_dim2 = dim1_or_dim2.value
    if dim1_or_dim2 == 'mean':
        return torch.load(pt_file)['mean_representations'][repr_layer]
    elif dim1_or_dim2 == 'per_tok':
        return torch.load(pt_file)['representations'][repr_layer]


def get_esm_embeddings(fasta_file, output_direc, dim1_or_dim2: EmbedDims, esm_model='esm2_t30_150M_UR50D'):
    if not os.path.exists(output_direc):
        os.makedirs(output_direc)
    os.system(
        rf'python C:\Users\jamel\PycharmProjects\Jamel\esm\scripts\extract.py {esm_model} {fasta_file}  {output_direc} --include {dim1_or_dim2.value}')



# def score_embedding_distance(parent_1_pt: str, parent_2_pt: str, chimera_pt: str, dim1_or_dim2: EmbedDims,
#                              repr_layer=30):
#     parent1_embed = pt_to_tensor(parent_1_pt, dim1_or_dim2, repr_layer)
#     parent2_embed = pt_to_tensor(parent_2_pt, dim1_or_dim2, repr_layer)
#     chi_embed = pt_to_tensor(chimera_pt, dim1_or_dim2, repr_layer)
#     return (abs(torch.sum(chi_embed - parent1_embed)) + abs(torch.sum(chi_embed - parent2_embed))).item()
def score_embedding_distance(input_fasta,):


def score_all_embeddings(pt_directory,  input_fasta: str, esm_model='esm2_t30_150M_UR50D'):
    def label_to_file(label):
        return os.path.join(pt_directory, f'{label}.pt')
    distance_scores = pd.DataFrame(columns=['chimera', 'score'])
    last_embed_layer = int(re.search(r't(\d*)_', esm_model).group(1))
    chimera_ids = create_dictionary_from_alignment(input_fasta).keys()
    distance_scores['chimera'] = chimera_ids
    for row, data in distance_scores.iterrows():
        parent1, parent2, _, _ =data['chimera'].split('-')
        distance_scores.loc[row, 'score'] = score_embedding_distance(
            label_to_file(parent1), label_to_file(parent2), label_to_file(data['chimera']), EmbedDims.Dim1,
            last_embed_layer)
    return distance_scores


def embedding_umap(pt_directory, parent_labels: tuple, input_fasta: str, esm_model='esm2_t30_150M_UR50D'):
    def label_to_file(label):
        return os.path.join(pt_directory, f'{label}.pt')

    parent1, parent2 = parent_labels
    distance_scores = pd.DataFrame(columns=['chimera', 'score'])
    chimera_ids = create_dictionary_from_alignment(input_fasta).keys()
    distance_scores['chimera'] = chimera_ids
    embeddings=[]
    for row, data in distance_scores.iterrows():
        embeddings.append(pt_to_tensor(label_to_file(data['chimera']), EmbedDims.Dim1))
        distance_scores.loc[row, 'score'] = score_embedding_distance(
            label_to_file(parent1), label_to_file(parent2), label_to_file(data['chimera']), EmbedDims.Dim1)

    #colorcode by sequence similarity?
    reducer = umap.UMAP(n_components=2)
    # rows in umap matrix are samples/proteins
    embedding_matrix = np.vstack(tuple(embeddings))
    embedding = reducer.fit_transform(embedding_matrix)
    plt.scatter(embedding[:, 0], embedding[:, 1], c=distance_scores['score'], cmap='viridis',marker='+')
    plt.scatter(embedding[:2,0], embedding[:2,1], color='red', edgecolor='black', s=150, marker='+', label='parents')
    plt.title('Riff CsChrim UMAP')
    plt.xlabel('Component 1')
    plt.ylabel('Component 2')
    plt.colorbar(label='Embedding Distance')
    plt.show()


class ESMArguments:
    parent_aln_file: str
    '''An fasta alignment file with your two parent proteins'''
    splice_size: int
    '''The size of the splice site that will move down the alignmentat some step size x and recombine the two parent sequences'''
    chimera_seq_output: str
    '''A single fasta file that contains all the chimeric sequences created by some splice and step size. (Including the parents)'''
    step_size: int = 1
    '''How much the splice site steps over after a recombination'''
    esm_output_directory: str
    '''Directory in create all esm outputs (.pt files) and the chimera_seq_output. If the directory doesn't exist it'll be created for you.'''
    embeddings_dimension: int = 1
    '''The ESM embeddings can be made per residue (2-dimensional) or just per sequence (1-dimensional)'''
    esm_model: str = 'esm2_t30_150M_UR50D'
    '''There are multiple different esm models with increasing amount of network layers and parameters. This 
    prediction models uses 150M parameter model for greatest performance (esm2_t48_15B_UR50D,esm2_t36_3B_UR50D,
    esm2_t33_650M_UR50D,esm2_t30_150M_UR50D,esm2_t12_35M_UR50D,esm2_t6_8M_UR50D)'''
    distance_score_csv: str
    '''Csv file that contains the embedding distances'''

    def __init__(self, arg_json):
        with open(arg_json, 'rb') as jfile:
            arg_dict = load(jfile)
        for key, value in arg_dict.items():
            setattr(self, key, value)

if __name__=='__main__':

    parser = argparse.ArgumentParser(
        description='Creating a fasta file of all potential chimeric proteins between two parents based on a sliding splice site')
    parser.add_argument('-in', '--inputjson', type=str, required=True,
                        help='alignment file in fasta style between 2 proteins')
    args = parser.parse_args()
    if args.inputjson:
        esm_args = ESMArguments(args.inputjson)
        parent_seq_dict = create_dictionary_from_alignment(esm_args.parent_aln_file)
        chi_combo_dict = create_chimera_combinations(parent_seq_dict, esm_args.splice_size,
                                                     scanner_rate=esm_args.step_size,
                                                     new_fasta_file=esm_args.chimera_seq_output)
        # parent_label = tuple(label.replace('-','') for label in parent_seq_dict.keys())
        if os.path.exists(esm_args.chimera_seq_output):
            # embedding_umap(esm_args.esm_output_directory,parent_label,esm_args.chimera_seq_output)
            get_esm_embeddings(esm_args.chimera_seq_output, esm_args.esm_output_directory, EmbedDims.Dim1)
            scores = score_all_embeddings(esm_args.esm_output_directory, esm_args.chimera_seq_output)
            # scores.to_csv(esm_args.distance_score_csv)
