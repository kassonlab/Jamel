import Immersion
from AccessiontoAlignment import create_dictionary_from_alignment, clustalw_to_fasta
# with open("/gpfs/gpfs0/scratch/jws6pq/Notebook/Overall/List_of_coronaviruses", 'r') as loc:
#     loc = loc.readlines()
#     label_pdb_list = []
#     chi_label_pdb_list = []
#     for line in loc:
#         file_stem = line.split()[-1]
#         label_pdb_list.append((file_stem, f'/gpfs/gpfs0/scratch/jws6pq/Notebook/PDB/3mer{file_stem}.pdb'))
#         chi_label_pdb_list.append((file_stem, f'/gpfs/gpfs0/scratch/jws6pq/Notebook/PDB/3merSARS2w{file_stem}S1.pdb'))
#         clustalw_to_fasta(f"/gpfs/gpfs0/scratch/jws6pq/Notebook/Alignment/{file_stem}onSARS2.aln",f"/gpfs/gpfs0/scratch/jws6pq/Notebook/Alignment/{file_stem}onSARS2.fasta")
# confidence_rank_matrix("/gpfs/gpfs0/scratch/jws6pq/Notebook/Fastas/SARS_S1_chimeras/CoronavirusMSA.aln",
#                        chi_label_pdb_list, '6vsb_B',
#                        '/gpfs/gpfs0/scratch/jws6pq/Notebook/Fastas/Immersion/rank_chi_plddt_matrix.csv',
#                        '/gpfs/gpfs0/scratch/jws6pq/Notebook/Fastas/Immersion/raw_chi_plddt_matrix.csv')
# confidence_rank_matrix("/gpfs/gpfs0/scratch/jws6pq/Notebook/Fastas/SARS_S1_chimeras/CoronavirusMSA.aln",
#                        label_pdb_list, '6vsb_B',
#                        '/gpfs/gpfs0/scratch/jws6pq/Notebook/Fastas/Immersion/rank_native_plddt_matrix.csv',
#                        '/gpfs/gpfs0/scratch/jws6pq/Notebook/Fastas/Immersion/raw_native_plddt_matrix.csv')

from numpy import zeros,savetxt


S1='QCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPGSASS'
alignment_file='/scratch/jws6pq/Notebook/Fixed_S1/Analysis_Folders/Fixed_S1/6vsb_MSA.aln'
ref_label='6VSB'
# native_rank='/scratch/jws6pq/Notebook/Fixed_S1/Analysis_Folders/Fixed_S1/fixed_native_ranked_matrix.npy'
# chimera_rank='/scratch/jws6pq/Notebook/Fixed_S1/Analysis_Folders/Fixed_S1/fixed_chi_ranked_matrix.npy'
# Immersion.confidence_rank_matrix(alignment_file,'/scratch/jws6pq/Notebook/Fixed_S1/Analysis_Folders/Fixed_S1/3mer6vsbw*S1.pdb',ref_label,'/scratch/jws6pq/Notebook/Fixed_S1/Analysis_Folders/Fixed_S1/fixed_chi_rank_matrix.csv',S1,chimera_rank,raw_matrix='/scratch/jws6pq/Notebook/Fixed_S1/Analysis_Folders/Fixed_S1/fixed_chi_raw_matrix.csv')
# Immersion.confidence_rank_matrix(alignment_file,'/scratch/jws6pq/Notebook/Fixed_S1/Analysis_Folders/Fixed_S1/3mer*.pdb',ref_label,'/scratch/jws6pq/Notebook/Fixed_S1/Analysis_Folders/Fixed_S1/fixed_native_rank_matrix.csv',S1,native_rank,'/scratch/jws6pq/Notebook/Fixed_S1/Analysis_Folders/Fixed_S1/fixed_native_raw_matrix.csv')
# TODO triple check
aln = create_dictionary_from_alignment(alignment_file)
aln_seq_length = len(aln[ref_label])
# correl_array=zeros((3,aln_seq_length))
# for residue in range(0,aln_seq_length):
#     correl_array[0,residue]=residue
#     correl_array[1, residue] = Immersion.residue_rank_correl(native_rank,chimera_rank,residue)[1]
#     correl_array[2, residue] = Immersion.residue_rank_correl(native_rank, chimera_rank, residue)[0]
# aln='/scratch/jws6pq/Notebook/Fixed_S1/Analysis_Folders/Fixed_S1/CoronavirusMSA.aln'
# comparison_matrix=zeros((len(aln),7),dtype=object)
# index=1239
# index_to_find=1765
# for indexes,label in enumerate(aln.keys()):
#     # protein_label=label.split()[-1]
#     protein_label=label
#     native=f'/scratch/jws6pq/Notebook/Fixed_S1/Analysis_Folders/Fixed_S1/3mer{protein_label}.pdb'
#     chi=f'/scratch/jws6pq/Notebook/Fixed_S1/Analysis_Folders/Fixed_S1/3mer6vsbw{protein_label}S1.pdb'
#     comparison=Immersion.contact_contingency(alignment_file,native,chi,'B',label,index,index_to_find,'6VSB',S1)
#     print(comparison)
#     comparison_matrix[indexes]=comparison
#
# savetxt(f'/scratch/jws6pq/Notebook/Fixed_S1/Analysis_Folders/Fixed_S1/1239_contingency.csv', comparison_matrix, fmt='%s/%s/%s/%s/%s/%s/%s')
# for 1379
comparison_matrix=zeros((len(aln),7),dtype=object)
index=1379
for indexes,label in enumerate(aln.keys()):
    # protein_label=label.split()[-1]
    protein_label=label
    native='/scratch/jws6pq/Notebook/Fixed_S1/Analysis_Folders/Fixed_S1/3mer{0}.pdb'
    chi='/scratch/jws6pq/Notebook/Fixed_S1/Analysis_Folders/Fixed_S1/3mer6vsbw{0}S1.pdb'
    comparison=Immersion.compare_aligned_contact_all_proteins(alignment_file,native,chi,'B',index,'/scratch/jws6pq/Notebook/Fixed_S1/Analysis_Folders/Fixed_S1/1379_contact_compare',S1)
    # print(comparison)
    # comparison_matrix[indexes]=comparison

# savetxt(f'/scratch/jws6pq/Notebook/Fixed_S1/Analysis_Folders/Fixed_S1/1239_contingency.csv', comparison_matrix, fmt='%s/%s/%s/%s/%s/%s/%s')
