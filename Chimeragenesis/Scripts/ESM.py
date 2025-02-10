import ast, regex, re, torch, umap, argparse, os
import difflib
from enum import Enum
from json import load
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle
from Bio import SeqIO
from AccessiontoAlignment import create_seq_records, fasta_creation
from ChimeraGenerator import sequence_splice, chimera_sequence_creation, \
    create_chimera_combinations
from prost import prost_embed_aln
from ankh_chimera import ankh_embed_aln
from Bio.Align import PairwiseAligner
from pandasgui import show


def word_match(word, choices, cutoff=0.5) -> str:
    return next(iter(difflib.get_close_matches(word, choices, n=1, cutoff=cutoff)), None)


class EmbedDims(Enum):
    Dim1 = 'mean'
    Dim2 = 'per_tok'

    def __str__(self):
        return self.name


class ESMModels(Enum):
    B15 = 'esm2_t48_15B_UR50D'
    B3 = 'esm2_t36_3B_UR50D'
    M605 = 'esm2_t33_650M_UR50D'
    M150 = 'esm2_t30_150M_UR50D'
    M35 = 'esm2_t12_35M_UR50D'
    M8 = 'esm2_t6_8M_UR50D'


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


r'C:\Users\jamel\PycharmProjects\Jamel\esm\scripts\extract.py'


def get_esm_embeddings(fasta_file, extract_py, output_direc, dim1_or_dim2: EmbedDims = EmbedDims.Dim2,
                       esm_model: ESMModels = ESMModels.M150):
    if not os.path.exists(output_direc):
        os.makedirs(output_direc)
    os.system(
        f'python {extract_py} {esm_model.value} {fasta_file}  {output_direc} --include {dim1_or_dim2.value}')
    output_direc = Path(output_direc)
    save_esm_embedding_directory(output_direc, output_direc.joinpath(output_direc.stem).with_suffix('.pkl'))


class SequenceDataframe(pd.DataFrame):
    def __init__(self, alignment_file, embed_dict_pkl: str = None):
        super().__init__(columns=['aln_sequence', 'sequence', 'description'])
        self.alignment_file = alignment_file
        self.create_dataframe_from_alignment()
        if embed_dict_pkl:
            self.EMBEDDINGS_DICT: dict[str,torch.Tensor] = torch.load(embed_dict_pkl, map_location=torch.device('cpu'))

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

    def get_description(self, chimera_label):
        if self.loc[chimera_label, 'description']:
            return ast.literal_eval(self.loc[chimera_label, 'description'])

    def add_protein(self, label, aln_seq, description):
        self.add_value(label, 'aln_sequence', aln_seq)
        self.add_value(label, 'sequence', aln_seq.replace('-', ''))
        self.add_value(label, 'description', description)

    def get_sequence(self, chimera_label):
        return self.loc[chimera_label, 'sequence']

    def get_aln(self, chimera_label):
        return self.loc[chimera_label, 'aln_sequence']

    def score_embedding_distance(self, chi_label, parent_labels, dist_func):
        distance = []
        for parent in parent_labels:
            distance.append(tensor_distance(self.EMBEDDINGS_DICT[chi_label], self.EMBEDDINGS_DICT[parent], dist_func))
        if dist_func != NormType.cosine:
            self.loc[chi_label, dist_func.name] = sum(distance)
        else:
            self.loc[chi_label, dist_func.name] = sum(distance)
        return sum(distance)

    def per_residue_distance_for_schema_set(self, chi_label):
        if self.EMBEDDINGS_DICT[chi_label].shape[0] != len(self.loc[chi_label, 'sequence']):
            return
        distance = 0

        for parent, aln_positions in residue_inheritance.items():
            for pos in aln_positions:
                distance += tensor_distance(self.EMBEDDINGS_DICT[chi_label][pos, :],
                                            self.EMBEDDINGS_DICT[parent][pos, :])
        return distance

    #TODO Use nested dictionaries or another data struct to show parent splice
    #TODO also need to account for mid-section splices
    def match_seq_to_parent(self, chi_label, parent, splice: slice):
        pattern = fr"{self.get_aln(chi_label)[splice]}"
        matches = regex.findall(pattern, self.get_aln(parent), overlapped=True)
        if len(matches) == 1:
            return slice(*re.search(pattern, self.get_aln(parent)).span())
        raise Exception("Couldn't find matching parent sequence")

    def score_per_res(self, chi_label, parent_labels: dict, dist_func):
        distance=0
        for parent, splice in parent_labels.items():
            chimera_slice = slice(*splice)
            parent_slice = self.match_seq_to_parent(chi_label, parent, chimera_slice)
            chi_tensor=self.EMBEDDINGS_DICT[chi_label][chimera_slice]
            parent_tensor=self.EMBEDDINGS_DICT[parent][parent_slice]
            distance+=tensor_distance(chi_tensor, parent_tensor, dist_func)
        if dist_func != NormType.cosine:
            self.loc[chi_label, dist_func.name] = distance
        else:
            self.loc[chi_label, dist_func.name] = distance
        return distance

    def score_all_per_res(self, dist_func: NormType):
        for label in self.EMBEDDINGS_DICT.keys():
            if len(parent_label := self.get_description(label)) == 2:
                self.score_per_res(label, parent_label, dist_func)

    def score_all_embeddings(self, dist_func: NormType):
        self.create_column(dist_func.name)
        for label in self.EMBEDDINGS_DICT.keys():
            if len(parent_label := self.get_description(label)) == 2:
                self.score_embedding_distance(label, parent_label, dist_func)

    def score_alnment_similarity(self):
        alignment = PairwiseAligner()
        for chimera, data in self.iterrows():
            for parent in self.get_description(chimera):
                self.loc[chimera, f'{parent}_similarity'] = alignment.score(self.get_aln(chimera),
                                                                            self.get_aln(parent)) / alignment.score(
                    self.get_aln(parent), self.get_aln(parent))


class EmbeddingAnalysis:
    provided_full_aln: bool
    distance_score_csv: str
    '''Csv file that contains the embedding distances'''
    aln_file: str
    splice_size: int
    chimera_seq_fasta: str
    step_size: int
    esm_output_directory: str
    embeddings_dimension: EmbedDims
    esm_model: ESMModels
    extract_py_file: str

    def __init__(self, arg_json):
        with open(arg_json, 'rb') as jfile:
            self.argument_dict = load(jfile)
        for key, value in self.argument_dict.items():
            setattr(self, key, value)
        self.embeddings_dimension = EmbedDims.Dim1 if self.embeddings_dimension == 1 else EmbedDims.Dim2
        closest_model = word_match(self.esm_model, [model.value for model in ESMModels], 0.7)
        if closest_model:
            self.esm_model = {model.value: model for model in ESMModels}[closest_model]
        else:
            print("Couldn't find ESM model")


def pt_to_tensor(pt_file, layer=None) -> torch.Tensor:
    tensor = torch.load(pt_file)[word_match('representation', torch.load(pt_file).keys())]
    if not layer:
        layer = sorted(list(tensor.keys()))[-1]
    return tensor[layer]


def add_row_of_zeroes(tensor: torch.Tensor):
    return torch.cat((torch.zeros(1, tensor.size(1)), tensor), dim=0)


def save_esm_embedding_directory(directory, new_pkl_file):
    embed_dict = {}
    for file in Path(directory).iterdir():
        embed_dict[file.stem] = pt_to_tensor(file)
    torch.save(embed_dict, new_pkl_file)


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
    schema_data = pd.read_csv(r"/scratch/jws6pq/Notebook/ESM/Schema_valdation/Chimeragenesis/schema_data.csv",
                              index_col='chimera_block_ID')
    combined_df = schema_data
    for func in dist_func:
        embeddings = SequenceDataframe(aln_file, embed_pkl_file)
        embeddings.score_all_embeddings(func)
        dist_type = func.name
        combined_df = pd.concat((combined_df, embeddings), axis=1)
        combined_df['exp_rank'] = combined_df['mKate_mean'].rank(ascending=False)
        combined_df[f'{dist_type}_rank'] = combined_df[dist_type].rank(ascending=func != NormType.cosine)
    combined_df = combined_df.loc[:, ~combined_df.columns.duplicated()]

    if new_distance_file:
        combined_df.to_csv(new_distance_file)
    return combined_df


def combine_w_randoms(embed_pkl_file, aln_file, dist_func: list[NormType], new_distance_file=''):
    randoms_data = pd.read_csv("..\Data\processed_randoms.csv", index_col='Protein')
    combined_df = randoms_data
    for func in dist_func:
        embeddings = SequenceDataframe(aln_file, embed_pkl_file)
        embeddings.score_all_per_res(func)
        dist_type = func.name
        combined_df = pd.concat((combined_df, embeddings), axis=1)
        combined_df[f'{dist_type}_rank'] = combined_df[dist_type].rank(ascending=func != NormType.cosine)
    combined_df = combined_df.loc[:, ~combined_df.columns.duplicated()]

    if new_distance_file:
        combined_df.to_csv(new_distance_file)
    return combined_df


if __name__ == '__main__':
    # TODO do per res with 3di sequence and cosine similarity
    # TODO turn esm embeddings into pkl and put in rivanna
    # randoms = SequenceDataframe('..\Data\PDTvPDK_chimera.aln', r'..\Data\PDTvPDK_prost.pkl')
    # randoms.score_per_res('PDTwPDK_152', {'PDK6':(0,51),'PDTH-2':(51,None)})
    # l1 = combine_w_schema(r"/scratch/jws6pq/Notebook/ESM/Schema_valdation/Chimeragenesis/15B_2d.pkl", '/scratch/jws6pq/Notebook/ESM/Schema_valdation/Chimeragenesis/labeled_schema_aln',
    #                       [NormType.cosine, NormType.manhattan, NormType.euclidean], '15B_2d.csv')
    mod = '_res'
    # print('esm')
    # l1 = combine_w_randoms(r"..\Data\150M_esm.pkl", '..\Data\PDTvPDK_chimera.aln',
    #                        [NormType.cosine, NormType.manhattan, NormType.euclidean], rf'..\Data\randoms_150M_esm{mod}.csv')
    print('prost')
    l1 = combine_w_randoms(r"..\Data\PDTvPDK_prost.pkl", '..\Data\PDTvPDK_chimera.aln',
                           [NormType.cosine, NormType.manhattan, NormType.euclidean], rf'..\Data\randoms_prost{mod}.csv')
    # print('ankh')
    # l1 = combine_w_randoms(r"..\Data\PDTvPDK_ankh.pkl", '..\Data\PDTvPDK_chimera.aln',
    #                        [NormType.cosine, NormType.manhattan, NormType.euclidean], rf'..\Data\randoms_ankh{mod}.csv')
    parser = argparse.ArgumentParser(
        description='Creating a fasta file of all potential chimeric proteins between two parents based on a sliding splice site')
    parser.add_argument('-in', '--inputjson', type=str, required=True,
                        help='alignment file in fasta style between 2 proteins')
    args = parser.parse_args()

    if args.inputjson:
        esm_args = EmbeddingAnalysis(args.inputjson)
        if esm_args.provided_full_aln:
            get_esm_embeddings(esm_args.aln_file, esm_args.extract_py_file, esm_args.esm_output_directory,
                               esm_args.embeddings_dimension, esm_args.esm_model)
            # SequenceDataframe(esm_args.aln_file)
        # if os.path.exists(esm_args.chimera_seq_fasta):
        #     esm_args.get_esm_embeddings()
        #     esm_args.create_embedding_container()
        #     embeddings = esm_args.score_all_embeddings()
        #     scores = score_all_embeddings(esm_args.esm_output_directory, esm_args.chimera_seq_output)
        #     scores.to_csv(esm_args.distance_score_csv)
