import os
import matplotlib.pyplot as plt
import torch
import umap
from Chimeragenesis.AccessiontoAlignment import create_dictionary_from_alignment, dictionary_to_fasta
from Chimeragenesis.ChimeraGenerator import sequence_splice, chimera_sequence_creation


# print(test_embeddings['representations'][30].shape)


two_sequence_aligment_file = "C:\Research\SARSvsMERS.txt"
sequence_dictionary = create_dictionary_from_alignment(two_sequence_aligment_file)


def create_chimera_combinations(two_seq_dict: dict, scanner_length, scanner_start=0, scanner_rate=1):
    """Takes two sequence dictionary and creates all possible chimeras with given splice length"""

    def scanning_chimera_generator(base_sequence, partner_seq):
        chimeric_seq_dict: dict[tuple, str] = {}
        splice_boundaries: list[tuple] = [(x, x + scanner_length) for x in
                                          range(scanner_start, len(base_sequence) - scanner_length, scanner_rate) if
                                          x + scanner_length < len(base_sequence)]
        for boundary in splice_boundaries:
            # This is the sequence to be spliced into by the partner seq, the residues between a given boundary have
            # been cut and replaced with a marker character '-' for the aligned sequence from the partner to replace it
            marked_base_sequence = sequence_splice(base_sequence, boundary)[1]
            partner_replacement = sequence_splice(partner_seq, boundary)[0]
            chimeric_seq = chimera_sequence_creation(partner_replacement, marked_base_sequence)
            chimeric_seq_dict[boundary] = chimeric_seq
        return chimeric_seq_dict

    seq1, seq2 = two_seq_dict.values()
    seq1_label, seq2_label = two_seq_dict.keys()
    final_chimera_dict: dict[str, dict[tuple, str]] = {
        seq1_label + 'w' + seq2_label: scanning_chimera_generator(seq1, seq2),
        seq2_label + 'w' + seq1_label: scanning_chimera_generator(seq2, seq1)}
    return final_chimera_dict


nested_chimeras_dict = create_chimera_combinations(sequence_dictionary, 600, scanner_rate=200)


def chimera_dict_to_multlple_chimera_fasta(parent_seq_dict,nested_chimera_dict: dict[str, dict[tuple, str]],new_fasta_file):
    sequence_dict_for_fasta = parent_seq_dict
    for combo, chimera_dict in nested_chimera_dict.items():
        for boundary, seq in chimera_dict.items():
            boundary = (str(x) for x in boundary)
            sequence_dict_for_fasta[f'{combo}_{"_".join(boundary)}'] = seq

    dictionary_to_fasta(sequence_dict_for_fasta, new_fasta_file)
def pt_to_tensor(pt_file,repr_layer=30):
    return torch.load(pt_file)['representations'][repr_layer]

# chimera_dict_to_multlple_chimera_fasta(sequence_dictionary,nested_chimeras_dict,'SARSvsMERS_600scan.fasta')
reducer = umap.UMAP(n_components=2)
esm_data=r'C:\Users\jamel\PycharmProjects\Jamel\esm\data'
mers_embeddings=pt_to_tensor(r'C:\Users\jamel\PycharmProjects\Jamel\esm\data\MERS_5X59.pt')
sars=pt_to_tensor(r'C:\Users\jamel\PycharmProjects\Jamel\esm\data\SARS2_7T9J.pt')
for file in os.listdir(esm_data):
    pt_file=os.path.join(esm_data,file)
    chimera=pt_to_tensor(pt_file)
    embedding = reducer.fit_transform(chimera)
    print(embedding)
    # Plot the UMAP representation
    plt.scatter(embedding[:, 0], embedding[:, 1])
    plt.title('UMAP Representation')
    plt.xlabel('Component 1')
    plt.ylabel('Component 2')
    plt.show()
    # score=abs(torch.sum(chimera-mers_embeddings))+abs(torch.sum(chimera-sars))
    # print(file,score.item())