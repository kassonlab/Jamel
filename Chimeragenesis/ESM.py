import ast, re, torch, umap, argparse, os
from enum import Enum
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle
from Bio import SeqIO
from AccessiontoAlignment import create_seq_records, fasta_creation
from ChimeraGenerator import sequence_splice, chimera_sequence_creation, \
    create_chimera_combinations


class EmbedDims(Enum):
    Dim1 = 'mean'
    Dim2 = 'per_tok'

    def __str__(self):
        return self.name


def manhattan_norm(tensor1: torch.Tensor, tensor2: torch.Tensor):
    tensor = tensor1 - tensor2
    return tensor.to(torch.float64).flatten().abs().sum().item()


def euclidean_norm(tensor1: torch.Tensor, tensor2: torch.Tensor):
    tensor = tensor1 - tensor2
    return torch.linalg.vector_norm(tensor.flatten()).item()


def cosine_similarity(tensor1: torch.Tensor, tensor2: torch.Tensor):
    tensor2 = tensor2.flatten()
    tensor1 = tensor1.flatten()
    return (torch.dot(tensor1, tensor2) / (
                torch.linalg.vector_norm(tensor1) * torch.linalg.vector_norm(tensor2))).item()


func_dict = {'euclidean_norm': euclidean_norm, 'manhattan_norm': manhattan_norm, 'cosine_similarity': cosine_similarity}


class NormType(Enum):
    euclidean = 'euclidean_norm'
    manhattan = 'manhattan_norm'
    cosine = 'cosine_similarity'

    def execute(self, tensor1, tensor2):
        return func_dict[self.value](tensor1, tensor2)


def pkl_a_file(obj, new_pkl_file):
    with open(new_pkl_file, "wb") as f:
        pickle.dump(obj, f)


def unpkl_a_file(pkl_file):
    with open(pkl_file, "rb") as f:
        return pickle.load(f)


def label_to_file(pt_directory, label):
    return os.path.join(pt_directory, f'{label}.pt')


def tensor_distance(tensor1, tensor2, dist_func=NormType.euclidean):
    return dist_func.execute(tensor1, tensor2)


class SequenceDataframe(pd.DataFrame):
    def __init__(self, alignment_file):
        super().__init__(columns=['aln_sequence', 'sequence', 'description'])
        self.alignment_file = alignment_file
        self.create_dataframe_from_alignment()

    def dataframe_to_aln(self, new_fasta_file):
        seq_records = [create_seq_records(str(seq_id), seq_info['aln_sequence'], str(seq_info['description'])) for
                       seq_id, seq_info in self.iterrows()]
        fasta_creation(new_fasta_file, seq_records)

    def create_dataframe_from_alignment(self):
        """Takes a fasta style alignment and makes a dictionary where the key is whatever signifier follows '>'
        and the value is the sequence with no spaces"""
        with open(self.alignment_file) as handle:
            for seq in SeqIO.parse(handle, "fasta"):
                self.loc[seq.id] = {'sequence': str(seq.seq).replace('-', ''), 'aln_sequence': str(seq.seq),
                                    'description': seq.description.replace(f'{seq.id} ', '')}
        return self

    def make_individual_fasta(self, output_direc, subunit_count):
        for protein, info in self.iterrows():
            fasta_creation(Path(output_direc).joinpath(str(protein)).with_suffix('.fasta'),
                           create_seq_records(protein, info['sequence'], info['description'], subunit_count))

    def create_column(self, column_name):
        self[column_name] = np.nan

    def add_value(self, label, column, value):
        self.loc[label, column] = value

    def get_parents(self, chimera_label):
        return ast.literal_eval(self.loc[chimera_label, 'description'])

    def add_protein(self, label, aln_seq, description):
        self.add_value(label, 'aln_sequence', aln_seq)
        self.add_value(label, 'sequence', aln_seq.replace('-', ''))
        self.add_value(label, 'description', description)

    def get_sequence(self, chimera_label):
        return self.loc[chimera_label, 'sequence']

    def get_aln(self, chimera_label):
        return self.loc[chimera_label, 'aln_sequence']


class EmbeddingAnalysis:
    distance_score_csv: str
    '''Csv file that contains the embedding distances'''
    EMBEDDINGS_DICT: dict[str:torch.Tensor]

    def __init__(self, embed_dict_pkl: str, aln_file):
        self.EMBEDDINGS_DICT: dict = torch.load(embed_dict_pkl, map_location=torch.device('cpu'))
        self.aln_df = SequenceDataframe(aln_file)

    def score_embedding_distance(self, chi_label, parent_labels, dist_func=None):
        distance = 0
        for parent in parent_labels:
            distance += tensor_distance(self.EMBEDDINGS_DICT[chi_label], self.EMBEDDINGS_DICT[parent], dist_func)
        if self.aln_df.get('distance') is not None and not self.aln_df.empty:
            self.aln_df.loc[chi_label, 'distance'] = distance
        return distance

    def per_residue_distance_for_val_set(self, chi_label):
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
        self.aln_df.create_column('distance')
        for label, chi_tensor in self.EMBEDDINGS_DICT.items():
            self.per_residue_distance_for_val_set(label, self.aln_df.get_parents(label))

    def score_all_embeddings(self, dist_func):
        self.aln_df.create_column('distance')
        for label, chi_tensor in self.EMBEDDINGS_DICT.items():
            if len(parent_label := self.aln_df.get_parents(label)) == 2:
                self.score_embedding_distance(label, parent_label, dist_func)


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


def save_embedding_directory(directory, new_pkl_file):
    embed_dict = {}
    for file in Path(directory).iterdir():
        embed_dict[file.stem] = torch.load(file)
    torch.save(embed_dict,new_pkl_file)

def embedding_umap(embed_pkl_file: str):
    reducer = umap.UMAP(n_components=2)
    embed_dict: dict[str, torch.Tensor] = torch.load(embed_pkl_file)
    embed_dict = {label: torch.sum(tensor, dim=0) for label, tensor in embed_dict.items()}
    # rows in umap matrix are samples/proteins
    embedding_matrix = np.vstack(tuple(embed_dict.values()))
    umap_vectors = reducer.fit_transform(embedding_matrix)
    plt.scatter(umap_vectors[:, 0], umap_vectors[:, 1], marker='+')
    plt.show()
    return plt.gcf()


def combine_w_schema(embed_pkl_file, aln_file, dist_func: list[NormType], new_distance_file=''):
    schema_data = pd.read_csv(r"/scratch/jws6pq/Notebook/ESM/Schema_valdation/Chimeragenesis/schema_data.csv", index_col='chimera_block_ID')
    combined_df = schema_data
    for func in dist_func:
        embeddings = EmbeddingAnalysis(embed_pkl_file, aln_file)
        embeddings.score_all_embeddings(func)
        dist_type = func.name
        combined_df = pd.concat((combined_df, embeddings.aln_df), axis=1)
        combined_df = combined_df.rename(columns={"distance": dist_type})
        combined_df['exp_rank'] = combined_df['mKate_mean'].rank(ascending=False)
        combined_df[f'{dist_type}_rank'] = combined_df[dist_type].rank(ascending=func!=NormType.cosine)
    combined_df = combined_df.loc[:, ~combined_df.columns.duplicated()]

    if new_distance_file:
        combined_df.to_csv(new_distance_file)
    return combined_df


if __name__ == '__main__':
    # save_embedding_directory('/scratch/jws6pq/Notebook/ESM/Schema_valdation/Chimeragenesis/15B_model/2d_schema_data/','/scratch/jws6pq/Notebook/ESM/Schema_valdation/Chimeragenesis/15B_2d.pkl')
    # TODO do per res with 3di sequence and cosine similarity
    # TODO turn esm embeddings into pkl and put in rivanna
    # l1 = combine_w_schema(r'full_schema_3di1.pkl', 'labeled_schema_aln',
    #                       [NormType.cosine,NormType.manhattan,NormType.euclidean],'full_schema_3di.csv')
    # l1 = combine_w_schema(r'large_ankh.pkl', 'labeled_schema_aln',
    #                       [NormType.cosine, NormType.manhattan, NormType.euclidean], 'large_ankh.csv')
    # l1 = combine_w_schema(r'prost.pkl', 'labeled_schema_aln',
    #                       [NormType.cosine, NormType.manhattan, NormType.euclidean], 'prost.csv')
    # l1 = combine_w_schema(r'alpha_full_3di.pkl', 'labeled_schema_aln',
    #                       [NormType.cosine, NormType.manhattan, NormType.euclidean], 'alpha_3di.csv')
    l1 = combine_w_schema(r"/scratch/jws6pq/Notebook/ESM/Schema_valdation/Chimeragenesis/15B_2d.pkl", '/scratch/jws6pq/Notebook/ESM/Schema_valdation/Chimeragenesis/labeled_schema_aln',
                          [NormType.cosine, NormType.manhattan, NormType.euclidean], '15B_2d.csv')
    # embedding_umap(r'alpha_full_3di.pkl').show()
    # parser = argparse.ArgumentParser(
    #     description='Creating a fasta file of all potential chimeric proteins between two parents based on a sliding splice site')
    # parser.add_argument('-in', '--inputjson', type=str, required=True,
    #                     help='alignment file in fasta style between 2 proteins')
    # args = parser.parse_args()
    # if args.inputjson:
    #     esm_args = EmbeddingAnalysis(args.inputjson)

    #     if os.path.exists(esm_args.chimera_seq_fasta):
    #         # esm_args.get_esm_embeddings()
    #         esm_args.create_embedding_container()
    #         embeddings = esm_args.score_all_embeddings()
    #         scores = score_all_embeddings(esm_args.esm_output_directory, esm_args.chimera_seq_output)
    #         scores.to_csv(esm_args.distance_score_csv)
