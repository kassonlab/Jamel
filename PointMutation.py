from ChimeraGenerator import chimeracls
from ColorCoding import create_dictionary_from_alignment
import argparse
aln='CoronavirusMSA.aln'
mutate_to='D'
seq_dict=create_dictionary_from_alignment(aln)
AMINO_ACIDS = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
