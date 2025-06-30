import ast,difflib,os,pickle,re,regex,torch
from datetime import date
from enum import Enum
from json import load
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
from AccessiontoAlignment import create_seq_records, fasta_creation, run_emboss_needle, calculate_sequence_identity


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


# A set of functions used to aggregate embedding distances in SequenceDataframe
evaluation_pairs = {min: max, np.mean: np.mean, max: min, sum: sum}


# make cpp function for when no gpu
def manhattan_norm(tensor1: torch.Tensor, tensor2: torch.Tensor):
    tensor = tensor1 - tensor2
    return tensor.flatten().to(torch.float64).abs().sum().item()


def euclidean_norm(tensor1: torch.Tensor, tensor2: torch.Tensor):
    tensor = tensor1 - tensor2
    return torch.linalg.vector_norm(tensor).item()


def cosine_similarity(tensor1: torch.Tensor, tensor2: torch.Tensor):
    tensor2 = tensor2.flatten()
    tensor1 = tensor1.flatten()
    return (torch.dot(tensor1, tensor2) / (
            torch.linalg.vector_norm(tensor1) * torch.linalg.vector_norm(tensor2))).item()


def dot_product(tensor1: torch.Tensor, tensor2: torch.Tensor):
    tensor2 = tensor2.flatten()
    tensor1 = tensor1.flatten()
    return torch.dot(tensor1, tensor2).item()


# def random_thing(tensor1: torch.Tensor, tensor2: torch.Tensor):
#     #TODO make a dot product that is truly
#
#     return c.max().item()


# Used in conjunction with NormType enum for distance function selection
func_dict = {'euclidean_norm': euclidean_norm, 'manhattan_norm': manhattan_norm, 'cosine_similarity': cosine_similarity,
              'dot_product': dot_product}


class NormType(Enum):
    euclidean = 'euclidean_norm'
    manhattan = 'manhattan_norm'
    cosine = 'cosine_similarity'
    random = 'random'
    dot_product = 'dot_product'

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


def get_esm_embeddings(fasta_file, extract_py, output_direc, dim1_or_dim2: EmbedDims = EmbedDims.Dim2,
                       esm_model: ESMModels = ESMModels.M150):
    if not os.path.exists(output_direc):
        os.makedirs(output_direc)
    os.system(
        f'python {extract_py} {esm_model.value} {fasta_file}  {output_direc} --include {dim1_or_dim2.value}')
    output_direc = Path(output_direc)
    save_esm_embedding_directory(output_direc, output_direc.joinpath(output_direc.stem).with_suffix('.pkl'))


class SequenceDataframe(pd.DataFrame):
    """The SequenceDataframe is a special dataframe used to organize and manipulate sequences from a fasta alignment.
    Especially for chimeras and in conjunction with llm embeddings. An empty can also be created to add freshly created sequences later.
    As well as spit out a fasta alignment."""

    def __init__(self, alignment_file=None, embed_dict_pkl: str = None, unconverted_df: pd.DataFrame = None):
        if isinstance(unconverted_df, pd.DataFrame):
            super().__init__(unconverted_df)
        else:
            super().__init__(columns=['aln_sequence', 'sequence', 'description'])
        device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
        self.emboss_format = '{0}w{1}.emb'
        self.alignment_file = alignment_file
        if alignment_file:
            self.create_dataframe_from_alignment()
        if embed_dict_pkl:
            self.EMBEDDINGS_DICT: dict[str, torch.Tensor] = torch.load(embed_dict_pkl, map_location=device)

    def dataframe_to_aln(self, new_fasta_file):
        seq_records = sum([create_seq_records(str(seq_id), seq_info['aln_sequence'], str(seq_info['description'])) for
                           seq_id, seq_info in self.iterrows()], [])
        fasta_creation(new_fasta_file, seq_records)

    def dataframe_to_multi_fa(self, new_fasta_file):
        seq_records = sum([create_seq_records(str(seq_id), seq_info['sequence'], str(seq_info['description'])) for
                           seq_id, seq_info in self.iterrows()], [])
        fasta_creation(new_fasta_file, seq_records)

    def create_dataframe_from_alignment(self):
        """Takes a fasta style alignment and makes a dictionary where the key is whatever signifier follows '>'
        and the value is the sequence with no spaces"""
        with open(self.alignment_file) as handle:
            for seq in SeqIO.parse(handle, "fasta"):
                self.loc[seq.id] = {'sequence': str(seq.seq).replace('-', ''), 'aln_sequence': str(seq.seq),
                                    'description': seq.description.replace(f'{seq.id} ', '')}
        return self

    def make_individual_fasta(self, output_direc, subunit_count=1):
        for protein, info in self.iterrows():
            fasta_creation(Path(output_direc).joinpath(str(protein)).with_suffix('.fa'),
                           create_seq_records(protein, info['sequence'], str(info['description']), subunit_count))

    def create_column(self, column_name):
        self[column_name] = np.nan

    def add_value(self, label, column, value):
        self.at[label, column] = value

    def get_description(self, chimera_label):
        # In this context, description should be used for labeling parents as for chimeras they appear in the alignment in dictionary style.
        # Current style is {parent:{parent_splice:corresponding_chimera_splice}}
        # you can use it for traditional descriptions, but you cannot then use the embedding scoring functions
        if self.loc[chimera_label, 'description']:
            return ast.literal_eval(self.loc[chimera_label, 'description'])
        return None

    def add_protein(self, label, aln_seq, description):
        self.loc[label] = {'sequence': str(aln_seq).replace('-', ''), 'aln_sequence': str(aln_seq),
                           'description': description}

    def get_sequence(self, chimera_label):
        return self.at[chimera_label, 'sequence']

    def get_aln(self, chimera_label):
        return self.at[chimera_label, 'aln_sequence']

    def score_embedding_distance(self, chi_label, parent_labels, dist_func, evaluation_func=min):
        distance = []
        for parent in parent_labels:
            distance.append(tensor_distance(self.EMBEDDINGS_DICT[chi_label], self.EMBEDDINGS_DICT[parent], dist_func))
        if dist_func == NormType.cosine or dist_func == NormType.dot_product:
            self.loc[chi_label, dist_func.name] = evaluation_pairs[evaluation_func](distance)
        else:
            self.loc[chi_label, dist_func.name] = evaluation_func(distance)
        return evaluation_func(distance)

    def score_all_embeddings(self, dist_func: NormType, evaluation_func=min):
        self.create_column(dist_func.name)
        for label in self.index:
            if len(parent_label := self.get_description(label)) == 2:
                self.score_embedding_distance(label, parent_label, dist_func, evaluation_func)

    def compare_confirmations(self, chi_label, parent_labels, dist_func, standards):
        possible_confirmations = {AA: [] for AA in "ARNDCEQGHILKMFPSTWYV"}
        distance = []

        for standard, tensor in standards.EMBEDDINGS_DICT.items():
            for residue, vector1 in zip(standards.get_aln(standard), tensor):
                if residue in possible_confirmations.keys():
                    possible_confirmations[residue].append(vector1)
        for residue, vector2 in zip(self.get_aln(chi_label), self.EMBEDDINGS_DICT[chi_label]):
            if residue in possible_confirmations.keys():
                res_confirmations = possible_confirmations[residue]
                if dist_func == NormType.cosine or dist_func == NormType.dot_product:
                    distance.append(max(tensor_distance(vector2, conf, dist_func) for conf in res_confirmations))
                else:
                    distance.append(min(tensor_distance(vector2, conf, dist_func) for conf in res_confirmations))
        distance = np.mean(distance)
        self.loc[chi_label, dist_func.name] = distance.detach()
        return distance

    def score_all_confs(self, dist_func: NormType, standards):
        for label in self.EMBEDDINGS_DICT.keys():
            if len(parent_label := self.get_description(label)) == 2:
                self.compare_confirmations(label, parent_label, dist_func, standards)

    def match_seq_to_parent(self, chi_label, parent, splice: slice):
        pattern = fr"{self.get_aln(chi_label)[splice]}"
        matches = regex.findall(pattern, self.get_aln(parent), overlapped=True)
        if len(matches) == 1:
            return slice(*re.search(pattern, self.get_aln(parent)).span())
        raise Exception("Couldn't find matching parent sequence")

    def score_per_res(self, chi_label, inheritance: dict[str, dict[tuple, tuple]], dist_func):
        distance = []
        for parent, splices in inheritance.items():
            for parent_splice, chi_splice in splices.items():
                # TODO need to ensure the shape matches with expectation, meaning embeddings are same length as sequences,
                # TODO although not sure how to check because some will have paddings,
                # TODO maybe i check for typical hidden dimension length??
                for parent_pos,chi_pos in zip(range(*parent_splice), range(*chi_splice)):
                    distance.append(tensor_distance(self.EMBEDDINGS_DICT[chi_label][chi_pos], self.EMBEDDINGS_DICT[parent][parent_pos], dist_func))
        self.loc[chi_label, dist_func.name] = np.mean(distance)
        return distance

    def score_all_per_res(self, dist_func: NormType):
        for label in self.EMBEDDINGS_DICT.keys():
            if len(inheritance := self.get_description(label)) == 2:
                self.score_per_res(label, inheritance, dist_func)

    def create_emboss_file(self, output_direc, needle_command='needle'):
        for chi in self.index:
            for parent in self.get_description(chi):
                emb_file = Path(output_direc).joinpath(self.emboss_format.format(chi, parent))
                emb_file = translate_windows_path(str(emb_file)) if str(emb_file).find('\\') != -1 else emb_file
                run_emboss_needle(emb_file, self.get_sequence(chi), self.get_sequence(parent), needle_command)

    def get_sequence_identity(self):
        for chi in self.index:
            identity = [calculate_sequence_identity(self.get_aln(chi), self.get_aln(parent)) for parent in
                        self.get_description(chi)]
            self.add_value(chi, 'identity', max(identity))

    def save_df(self,file_name,metadata:dict=None):
        date_={'date':date.today().strftime('%Y-%m-%d')}
        self.loc['metadata', 'sequence'] = ""
        if metadata:
            date_.update(metadata) if metadata else date_
        self.at['metadata', 'sequence'] = [date_]
        self.to_csv(file_name,index_label='label')


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


# def embedding_umap_2d(embed_pkl_file: str):
#     reducer = umap.UMAP(n_components=2)
#     embed_dict: dict[str, torch.Tensor] = torch.load(embed_pkl_file, map_location=torch.device('cpu'))
#     embed_dict = {label: torch.sum(tensor, dim=0) for label, tensor in embed_dict.items()}
#     # rows in umap matrix are samples/proteins
#     embedding_matrix = np.vstack(tuple(embed_dict.values()))
#     umap_vectors = reducer.fit_transform(embedding_matrix)
#     plt.scatter(umap_vectors[:, 0], umap_vectors[:, 1], marker='+')
#     for i in range(len(umap_vectors[:, 0])):
#         plt.text(umap_vectors[i, 0], umap_vectors[i, 1], f"{list(embed_dict.keys())[i]}", fontsize=10, ha="right",
#                  va="bottom")
#     plt.show()
#     return plt.gcf()


def embedding_umap_3d(embed_pkl_file: str):
    embed_dict: dict[str, torch.Tensor] = torch.load(embed_pkl_file, map_location=torch.device('cpu'))
    embed_dict = {label: torch.sum(tensor, dim=0) for label, tensor in embed_dict.items()}
    embedding_matrix = np.vstack(tuple(embed_dict.values()))
    umap_3d = umap.UMAP(n_components=3, random_state=42)
    X_umap = umap_3d.fit_transform(embedding_matrix)

    # Create 3D scatter plot
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')

    # Scatter plot with colors based on labels
    sc = ax.scatter(X_umap[:, 0], X_umap[:, 1], X_umap[:, 2], cmap='Spectral')
    for i in range(len(X_umap[:, 0])):
        ax.text(X_umap[i, 0], X_umap[i, 1], z=X_umap[i, 2], s=f"{list(embed_dict.keys())[i]}", fontsize=10, ha="right",
                va="bottom")
    # Add colorbar
    cbar = plt.colorbar(sc)
    cbar.set_label("Digit Label")

    # Labels and title
    ax.set_title("3D UMAP Projection of Digits Dataset")
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    # ax.set_zlabel("UMAP 3")
    #
    # plt.show()
    # return plt.gcf()


def combine_w_schema(embed_pkl_file, aln_file, dist_func: list[NormType], evaluation_func, new_distance_file=''):
    schema_data = pd.read_csv(r"/Chimeragenesis/Data/schema_csvs/schema_data.csv",
                              index_col='chimera_block_ID')

    combined_df = schema_data
    for func in dist_func:
        embeddings = SequenceDataframe(aln_file, embed_pkl_file)
        embeddings.score_all_embeddings(func, evaluation_func)
        embeddings.get_sequence_identity()
        dist_type = func.name
        combined_df = pd.concat((combined_df, embeddings), axis=1)
        combined_df['exp_rank'] = combined_df['mKate_mean'].rank(ascending=False)
        combined_df[f'{dist_type}_rank'] = combined_df[dist_type].rank(
            ascending=func != NormType.cosine and func != NormType.dot_product)
    combined_df = combined_df.loc[:, ~combined_df.columns.duplicated()]

    if new_distance_file:
        combined_df.to_csv(new_distance_file)
    return combined_df


def combine_w_randoms(embed_pkl_file, aln_file, dist_func: list[NormType], evaluation_func, new_distance_file=''):
    randoms_data = pd.read_csv(r"..\Data\Randoms\processed_randoms.csv", index_col='Protein')
    combined_df = randoms_data
    for func in dist_func:
        print(func)
        embeddings = SequenceDataframe(aln_file, embed_pkl_file)
        embeddings.score_all_confs(func, SequenceDataframe(
            r'C:\Users\jamel\PycharmProjects\Jamel\Chimeragenesis\Data\pdth_siblings.aln',
            r'C:\Users\jamel\PycharmProjects\Jamel\Chimeragenesis\llm_output\pdth_siblings.pkl'))
        dist_type = func.name
        combined_df = pd.concat((combined_df, embeddings), axis=1)
        combined_df[f'{dist_type}_rank'] = combined_df[dist_type].rank(
            ascending=func != NormType.cosine and func != NormType.dot_product)
    combined_df = combined_df.loc[:, ~combined_df.columns.duplicated()]
    if new_distance_file:
        combined_df.to_csv(new_distance_file)
    return combined_df


def translate_windows_path(windows_path: str):
    drive = re.match(r'(\w):', windows_path)
    windows_path = windows_path.replace('\\', '/')
    if drive:
        windows_path = windows_path.replace(drive.group(0), '/mnt/' + drive.group(1).lower())
    return windows_path

if __name__ == '__main__':
#     # TODO do per res with 3di sequence and cosine similarity
#     # TODO turn esm embeddings into pkl and put in rivanna
#     save_esm_embedding_directory(r'C:\Users\jamel\PycharmProjects\Jamel\esm\3B_model\2d_schema_data',
#                                  '../llm_output/schema_3Besm.pkl')
#     # look at plot of sequence identity versus expression and see where there are ones that express high but are disimilar
#     schema_data = pd.read_csv(r"C:\Users\jamel\PycharmProjects\Jamel\Chimeragenesis\Data\schema_data.csv",
#                               index_col='chimera_block_ID')
#     notag=SequenceDataframe(r'C:\Users\jamel\PycharmProjects\Jamel\Chimeragenesis\Data\schema_no_tag.aln')
#     notag.get_sequence_identity()
#     notag=pd.concat((notag, schema_data), axis=1)
#     notag=notag.dropna()
#     notag=notag.dropna(axis=1)
#     plt.scatter(notag['identity'].astype(float),notag['mKate_mean'].astype(float))
#     for label,data in notag.iterrows():
#         plt.text(float(data['identity']),float(data['mKate_mean']), str(label), fontsize=10, ha="right",
#                  va="bottom")
#     plt.xlabel('Sequence Identity')
#     plt.ylabel('Expression (mKate)')
#     plt.title('Contiguous Chimera Expression vs. Identity')
#     plt.show()
    mod = 'mean'
    # l1=combine_w_schema(r"C:\Users\jamel\PycharmProjects\Jamel\Chimeragenesis\llm_output\schema_notag_X_v2.pkl", r'C:\Users\jamel\PycharmProjects\Jamel\Chimeragenesis\Data\schema_no_tag.aln',
    #                       [NormType.cosine, NormType.manhattan, NormType.euclidean,NormType.dot_product],
    #                       np.mean,f'../Data/notag_perres_dot.csv')
#     # combine_w_randoms(r"..\Data\150M_esm.pkl", '..\Data\PDTvPDK_chimera.aln',
#     #                        [NormType.cosine, NormType.manhattan, NormType.euclidean], rf'..\Data\randoms_150M_esm{mod}.csv')
#     # print('prost')
#     # combine_w_randoms(r"C:\Users\jamel\PycharmProjects\Jamel\Chimeragenesis\llm_output\randoms_prost_X.pkl",
#     #                   '..\Data\Randoms\PDTvPDK_chimera.aln',
#     #                   [NormType.cosine, NormType.manhattan, NormType.euclidean, NormType.dot_product],
#     #                   sum, rf'..\Data\randoms_prost_{mod}.csv')
    print('ankh')
    # combine_w_schema(r"C:\Users\jamel\PycharmProjects\Jamel\Chimeragenesis\llm_output\ankh_tensors.pkl", r'C:\Users\jamel\PycharmProjects\Jamel\Chimeragenesis\Data\schema_no_tag.aln',
    #                      [NormType.cosine, NormType.manhattan, NormType.euclidean, NormType.dot_product], np.mean,rf'..\Data\ankh_w_dot.csv')
    combine_w_schema(r"C:\Users\jamel\PycharmProjects\Jamel\Chimeragenesis\llm_output\schema_3Besm.pkl",
                     r'/Chimeragenesis/Data/schema_csvs/schema_no_tag.aln',
                     [NormType.cosine, NormType.manhattan, NormType.euclidean, NormType.dot_product], np.mean,rf'..\Data\esm_3B_w_dot.csv')
#     # parser = argparse.ArgumentParser(
#     #     description='Creating a fasta file of all potential chimeric proteins between two parents based on a sliding splice site')
#     # parser.add_argument('-in', '--inputjson', type=str, required=True,
#     #                     help='alignment file in fasta style between 2 proteins')
#     # args = parser.parse_args()
#     #
#     # if args.inputjson:
#     #     esm_args = EmbeddingAnalysis(args.inputjson)
#     #     if esm_args.provided_full_aln:
#     #         get_esm_embeddings(esm_args.aln_file, esm_args.extract_py_file, esm_args.esm_output_directory,
#     #                            esm_args.embeddings_dimension, esm_args.esm_model)
#     # SequenceDataframe(esm_args.aln_file)
#     # if os.path.exists(esm_args.chimera_seq_fasta):
#     #     esm_args.get_esm_embeddings()
#     #     esm_args.create_embedding_container()
#     #     embeddings = esm_args.score_all_embeddings()
#     #     scores = score_all_embeddings(esm_args.esm_output_directory, esm_args.chimera_seq_output)
#     #     scores.to_csv(esm_args.distance_score_csv)
