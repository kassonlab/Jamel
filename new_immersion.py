from Chimeragenesis.Analysis import get_info_from_plddt_file
from os import path,listdir
from Chimeragenesis.AccessiontoAlignment import create_dictionary_from_alignment,extract_seq_from_fasta
import pandas as pd
def raw_confidence_alignment(native_naming:str,chimera_naming:str,alignment_file:str,plddt_folder:str,recombination_sequence_fasta,):
    plddt_files=tuple(path.join(plddt_folder,file) for file in listdir(plddt_folder))
    native_scores={}
    chimera_scores={}
    aln_dict=create_dictionary_from_alignment(alignment_file)
    ref_recombined_site=extract_seq_from_fasta(recombination_sequence_fasta)

    for name,alnment in aln_dict.items():
        native_name=native_naming.replace('*',name)
        native_file=next((file for file in plddt_files if native_name in file),None)
        native_score=iter(get_info_from_plddt_file(native_file)[1])
        native_scores[name]=tuple(next(native_score,None) if acid.isalpha() else acid for acid in alnment)

        chimera_name = chimera_naming.replace('*', name)
        chimera_file = next((file for file in plddt_files if chimera_name in file), None)
        #Some iterators are getting exhausted before they should, got to find out why, Also need to truncate to S1
        chimera_score = iter(get_info_from_plddt_file(chimera_file)[1])
        chimera_scores[name] = tuple(next(chimera_score,None) if acid.isalpha() else acid for acid in alnment)

    print(pd.DataFrame(chimera_scores))
def replace_residue_w_score(recombined_sequence:str,):

if __name__=='__main__':
    raw_confidence_alignment('3mer*','3mer6vsbw*S1','6vsb_MSA.aln',r"C:\Research\Plddt",'Full_6vsb_S1.fasta')