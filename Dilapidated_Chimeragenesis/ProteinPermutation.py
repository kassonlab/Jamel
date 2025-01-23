import random


def permutate_protein(protein_seq:str):
    seq_list=list(protein_seq)
    random.shuffle(seq_list)
    return ''.join(seq_list)


def create_random_protein(protein_len:int):
    amino_acids=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L','M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    return ''.join(random.choices(amino_acids, k=protein_len))

