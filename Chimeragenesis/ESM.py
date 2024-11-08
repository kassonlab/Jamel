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


def label_to_file(pt_directory, label):
    return os.path.join(pt_directory, f'{label}.pt')


class ESMArguments:
    parent_aln_file: str
    '''An fasta alignment file with your two parent proteins'''
    splice_size: int
    '''The size of the splice site that will move down the alignment at some step size x and recombine the two parent sequences'''
    chimera_seq_fasta: str
    '''A single fasta file that contains all the chimeric sequences created by some splice and step size. (Including the parents)'''
    step_size: int = 1
    '''How much the splice site steps over after a recombination'''
    esm_output_directory: str
    '''Directory in create all esm outputs (.pt files) and the chimera_seq_output. If the directory doesn't exist it'll be created for you.'''
    embeddings_dimension = 1
    '''The ESM embeddings can be made per residue (2-dimensional) or just per sequence (1-dimensional)'''
    esm_model: str = 'esm2_t30_150M_UR50D'
    '''There are multiple different esm models with increasing amount of network layers and parameters. This 
    prediction models uses 150M parameter model for greatest performance (esm2_t48_15B_UR50D,esm2_t36_3B_UR50D,
    esm2_t33_650M_UR50D,esm2_t30_150M_UR50D,esm2_t12_35M_UR50D,esm2_t6_8M_UR50D)'''
    distance_score_csv: str
    '''Csv file that contains the embedding distances'''
    extract_py_file: str
    '''Path to the python script file from esm scripts called extract.py'''

    class EmbeddingLabel:
        def __init__(self, pt_file, dimensions: EmbedDims):
            self.pt_file = pt_file
            self.dimensions = dimensions
            self.label = Path(pt_file).stem
            self.embedding = pt_to_tensor(pt_file, dimensions)
            self.embedding = self.add_beg_of_seq(self.embedding)
            self.distance = None
            self.is_parent = False
            if '-' in self.label:
                self.base_protein, self.partner_protein, *self.splice_site = self.label.split('-')
                self.splice_site = tuple(int(x) for x in self.splice_site)
            else:
                self.is_parent = True

        def add_info_to_pt(self):
            pt = torch.load(self.pt_file)
            pt['base_protein'] = self.base_protein
            pt['partner_protein'] = self.partner_protein
            pt['splice_site'] = self.splice_site
            torch.save(pt, self.pt_file)

        def tensor_distance(self, tensor1, tensor2):
            return abs(torch.norm(tensor2,tensor1))

        def score_embedding_distance(self):
            if self.is_parent:
                return 0
            pt_direc = Path(self.pt_file).parent
            boundaries = tuple({0, self.splice_site[0], self.splice_site[1], self.embedding.shape[0]})
            boundaries = sorted(boundaries)
            boundaries = tuple((boundaries[index], boundaries[index + 1]) for index in range(len(boundaries) - 1))
            self.distance = 0

            base_pt = pt_to_tensor(label_to_file(pt_direc, self.base_protein), self.dimensions)
            partner_pt = pt_to_tensor(label_to_file(pt_direc, self.partner_protein), self.dimensions)
            base_pt = add_beg_of_seq(base_pt)
            partner_pt = add_beg_of_seq(partner_pt)
            for boundary in boundaries:
                if boundary == self.splice_site:
                    self.distance += self.tensor_distance(partner_pt[slice(*boundary)],
                                                          self.embedding[slice(*boundary)])
                else:
                    self.distance += self.tensor_distance(base_pt[slice(*boundary)], self.embedding[slice(*boundary)])
            self.distance = self.distance.item()
            return self.distance

    def __init__(self, arg_json):
        with open(arg_json, 'rb') as jfile:
            arg_dict = load(jfile)
        for key, value in arg_dict.items():
            setattr(self, key, value)
        self.embedding_container: list[ESMArguments.EmbeddingLabel] = []
        self.embeddings_dimension = EmbedDims.Dim1 if self.embeddings_dimension == 1 else EmbedDims.Dim2

    def create_embedding_container(self):
        for pt in os.listdir(self.esm_output_directory):
            if pt.endswith('pt'):
                self.embedding_container.append(
                    ESMArguments.EmbeddingLabel(os.path.join(self.esm_output_directory, pt), self.embeddings_dimension))

    def score_all_embeddings(self):
        return {embedding.label: embedding.score_embedding_distance() for embedding in self.embedding_container}

    def get_esm_embeddings(self):
        if not os.path.exists(self.esm_output_directory):
            os.makedirs(self.esm_output_directory)
        os.system(
            rf'python {self.extract_py_file} {self.esm_model} {self.chimera_seq_fasta}  {self.esm_output_directory} --include {self.embeddings_dimension.value}')
        self.create_embedding_container()


def pt_to_tensor(pt_file, dim1_or_dim2: EmbedDims = EmbedDims.Dim1, model='esm2_t30_150M_UR50D') -> torch.Tensor:
    dim1_or_dim2 = dim1_or_dim2.value
    repr_layer=int(re.search(r't(\d*)_',model).group(1))
    if dim1_or_dim2 == 'mean':
        tensor = torch.load(pt_file)['mean_representations'][repr_layer]
    else:
        tensor = torch.load(pt_file)['representations'][repr_layer]
    return tensor


def add_beg_of_seq(tensor: torch.Tensor):
    return torch.cat((torch.tensor([0]), tensor))


def embedding_umap(pt_directory, parent_labels: tuple, input_fasta: str, esm_model='esm2_t30_150M_UR50D'):
    parent1, parent2 = parent_labels
    distance_scores = pd.DataFrame(columns=['chimera', 'score'])
    chimera_ids = create_dictionary_from_alignment(input_fasta).keys()
    distance_scores['chimera'] = chimera_ids
    embeddings = []
    for row, data in distance_scores.iterrows():
        embeddings.append(pt_to_tensor(label_to_file(data['chimera']), EmbedDims.Dim1))
        distance_scores.loc[row, 'score'] = score_embedding_distance(
            label_to_file(parent1), label_to_file(parent2), label_to_file(data['chimera']), EmbedDims.Dim1)

    #colorcode by sequence similarity?
    reducer = umap.UMAP(n_components=2)
    # rows in umap matrix are samples/proteins
    embedding_matrix = np.vstack(tuple(embeddings))
    embedding = reducer.fit_transform(embedding_matrix)
    plt.scatter(embedding[:, 0], embedding[:, 1], c=distance_scores['score'], cmap='viridis', marker='+')
    plt.scatter(embedding[:2, 0], embedding[:2, 1], color='red', edgecolor='black', s=150, marker='+', label='parents')
    plt.title('Riff CsChrim UMAP')
    plt.xlabel('Component 1')
    plt.ylabel('Component 2')
    plt.colorbar(label='Embedding Distance')
    plt.show()


if __name__ == '__main__':

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
                                                     new_fasta_file=esm_args.chimera_seq_fasta)
        if os.path.exists(esm_args.chimera_seq_fasta):
            # esm_args.get_esm_embeddings()
            esm_args.create_embedding_container()
            embeddings = esm_args.score_all_embeddings()
            # scores = score_all_embeddings(esm_args.esm_output_directory, esm_args.chimera_seq_output)
            # scores.to_csv(esm_args.distance_score_csv)
