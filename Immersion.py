import math
from Bio import Seq
import Analysis
from Chimeragenesis.AccessiontoAlignment import get_alignment_indexing, no_gap_sequence_from_alignment, alignment_finder,create_dictionary_from_alignment
from numpy import zeros, savetxt, save,load
from scipy.stats import rankdata,spearmanr
import MDContacts
import ContactMap
from pathlib import Path

'''This is very specific code for dilapidated and wrong sectioning of S1 and conversion of individual alignments to 
map on a MSA'''

Sixvsb_S1 = 'QCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPGSASS'




def alignment_to_confidence(aln_file, comparison_label, ref_label, pdb_file, soi,include_AA=False):
    """Uses an MSA, the alignment label of your protein of interest, the label of your reference protein, the pdb of
    your comparison protein, the sequence of interest (soi) from the reference protein where the chimeric recombination
    occured, and returns a list the length of the alignment sequences, where all amino acids inside the soi are replaced
    with their corresponding confidence score according to their pdb.
    This only works with homologous proteins for now"""
    sequence_dict = create_dictionary_from_alignment(aln_file)
    comparison_alignment = sequence_dict[comparison_label]
    comparison_alignment_indexing = get_alignment_indexing(comparison_alignment)
    comparison_sequence = no_gap_sequence_from_alignment(comparison_alignment)
    comparison_soi = alignment_finder(aln_file, soi, comparison_label, ref_label)[0]
    # = ref_sequence.replace(seq_of_interest, homologous_splice)
    # Turns the alignment sequence into a mutable object (list) so that the corresponding residues can be swapped
    listed_alignment_characters = list(comparison_alignment)
    # Finds the no gap sequence indexes for the sequence of interest (soi)
    no_gap_start = comparison_sequence.find(comparison_soi)
    no_gap_end = comparison_sequence.find(comparison_soi) + len(comparison_soi)
    # Splices the indexing list into a list of indexes for residues that need to be replaced with the corresponding confidence score
    found_alignment_indexing = comparison_alignment_indexing[no_gap_start:no_gap_end]
    plddt = Analysis.get_plddt_dict_from_pdb(pdb_file)
    for sequence, plddt in plddt.items():
        start = sequence.find(comparison_soi)
        end = start + len(comparison_soi)
        plddt_chunk = plddt[start:end]
        for index, position in enumerate(found_alignment_indexing):
            listed_alignment_characters[position] = plddt_chunk[index]
    # If include boolean is true, all alignment positions outside the sequence of interest will remain untouched and
    # leave the respective amino acid at that position. Otherwise, all non-numbers are replaced with '-'
    if include_AA:
        return listed_alignment_characters
    confidence_scores = [x if not str(x).isalpha() and x!='-' else '-' for x in listed_alignment_characters]
    return confidence_scores
print(alignment_to_confidence('6vsb_MSA.aln','BAT2008','6VSB','3merBAT2008.pdb',Sixvsb_S1))

# turn this into smaller so you can put soi for each indivdual
def confidence_rank_matrix(alignment_file, pdb_naming, reference_label, rank_matrix_file, sequence_of_interest,
                           npy_file='',raw_matrix=''):
    #label columns with alignment index
    aln_dict = create_dictionary_from_alignment(alignment_file)
    aln_seq_length = len(aln_dict[reference_label])
    aln_length = len(aln_dict)
    plddt_matrix = zeros((aln_length + 1, aln_seq_length + 1), dtype=object)
    plddt_matrix[0,1:]=range(0,aln_seq_length)
    for index, label in enumerate(aln_dict.keys()):
        plddt_matrix[index+1, 0] = label
        plddt_matrix[index+1, 1:] = alignment_to_confidence(alignment_file, label, reference_label,
                                                          pdb_naming.replace('*', label), sequence_of_interest)
    if raw_matrix:
        savetxt(raw_matrix, plddt_matrix, fmt='%s')
    for index, column in enumerate(plddt_matrix[0, 1:]):
        column = plddt_matrix[1:, 1 + index]
        rank_column = rankdata([x for x in column if x != '-'], 'ordinal')
        column_indexing = [ind for ind, x in enumerate(column) if x != '-']
        for rank_index, correct_index in enumerate(column_indexing):
            column[correct_index] = rank_column[rank_index]
        plddt_matrix[1:, 1 + index] = column
    if npy_file:
        save(npy_file,plddt_matrix)


    savetxt(rank_matrix_file, plddt_matrix, fmt='%s')


def Truncating_sars_dna_sequence(aligned_full_length, aligned_trunc, unaligned_dna):
    codons = [unaligned_dna[index:index + 3] for index in range(0, len(unaligned_dna) - 3, 3)]
    Seq.Seq(unaligned_dna).translate()
    alignment_index_trunc = [ind for ind, x in enumerate(aligned_trunc) if x.isalpha()]
    trunc_dna = ''
    for index in alignment_index_trunc:
        if aligned_full_length.isalpha() and index != 985 and index != 986:
            trunc_dna += codons[index]
        # assigned matching proline codons based on wobble position of previous alignment
        elif index == 985:
            trunc_dna += 'CCA'
        elif index == 986:
            trunc_dna += 'CCT'
    return trunc_dna


def contact_contingency(alignment_file, native_pdb, chimera_pdb, chain_id, label, alignment_index,
                        contact_to_find_aln_index,ref_label,sequence_of_interest):
    try:
        residue = ContactMap.get_residue_at_native_position(alignment_file, label, alignment_index)
        native_index = ContactMap.correct_alignment_for_residue_position(alignment_file, label, alignment_index)
    except:
        print('doesnt exist')
        return [label, '-', '-', '-', '-', '-', '-']
    # TODO make a saveable inter
    native_contacts = MDContacts.get_intra_residue_contact_list(native_pdb,chain_id,Path(native_pdb).stem+'_matrix')[native_index]
    native_contacts += MDContacts.inter_residue_contact_list_md_analysis(native_pdb, chain_id)[native_index]
    native_contacts = [int(str(contacts).split(':')[-1]) for contacts in native_contacts]
    try:
        native_index_to_find = ContactMap.correct_alignment_for_residue_position(alignment_file, label,
                                                                                 contact_to_find_aln_index)
        native_yes = native_index_to_find in native_contacts
        native_no = native_index_to_find not in native_contacts
        native_seq = ContactMap.get_sequence_from_pdb(native_pdb, chain_id)
        residue_to_find = native_seq[
            ContactMap.correct_alignment_for_residue_position(alignment_file, label, contact_to_find_aln_index)]
    except:
        native_no = True
        native_yes = False
        residue_to_find = '-'

    chi_index = ContactMap.correct_alignment_for_chimera_index(alignment_file, ref_label, alignment_index, label,sequence_of_interest)[0]
    residue_ids = f'{residue}{native_index + 1},{residue}{chi_index + 1}'
    chimera_contacts = \
    MDContacts.get_intra_residue_contact_list(chimera_pdb, chain_id, Path(chimera_pdb).stem+'_matrix')[
        chi_index]
    chimera_contacts += ContactMap.get_inter_protein_contacts(chimera_pdb, chain_id)[chi_index]
    chimera_contacts = [int(str(contacts).split(':')[-1]) for contacts in chimera_contacts]
    chi_index_to_find = ContactMap.correct_alignment_for_chimera_index(alignment_file, ref_label, contact_to_find_aln_index, label,sequence_of_interest)[0]
    chi_yes = chi_index_to_find in chimera_contacts
    chi_no = chi_index_to_find not in chimera_contacts

    return [label, residue_ids, int(native_yes), int(native_no), int(chi_yes), int(chi_no),
            residue_to_find]

def create_rank_difference(native_rank_npy_file,chimera_rank_npy_file,aln_index,rank_difference_npy,rank_difference_file=''):
    native_matrix = load(native_rank_npy_file, allow_pickle=True)
    chimera_matrix = load(chimera_rank_npy_file, allow_pickle=True)
    difference_table=zeros((native_matrix.shape[0]-1,2),dtype=object)
    for index,name in enumerate(native_matrix[1:,0]):
        difference_table[index,0]=name
        if chimera_matrix[index+1,aln_index+1]=='-' or chimera_matrix[index+1,aln_index+1]=='-':
            difference_table[index, 1] = '-'
        else:
            difference_table[index,1]=chimera_matrix[index+1,aln_index+1]-native_matrix[index+1,aln_index+1]
    save(rank_difference_npy,difference_table)
    if rank_difference_file:
        savetxt(rank_difference_file,difference_table,fmt='%s')


def residue_rank_correl(native_rank_npy_file,chimera_rank_npy_file,residue_index):
    native_matrix=load(native_rank_npy_file,allow_pickle=True)
    chimera_matrix = load(chimera_rank_npy_file,allow_pickle=True)
    native_no_gap_column=[x for x in native_matrix[1:,residue_index+1] if x!='-']
    chimera_no_gap_column=[x for x in chimera_matrix[1:,residue_index+1] if x!='-']
    residue_saturation=len(native_no_gap_column)
    return spearmanr(native_no_gap_column,chimera_no_gap_column)[0],residue_saturation


# with open("/gpfs/gpfs0/scratch/jws6pq/Notebook/Overall/List_of_coronaviruses", 'r') as loc:
#     loc = loc.readlines()
# aln='/gpfs/gpfs0/scratch/jws6pq/Notebook/PDB/CoronavirusMSA.aln'
# index=1792
# soi='SSNF'
# comparison_matrix=zeros((len(loc),3),dtype=object)
# for indexes,label in enumerate(loc):
#     protein_label=label.split()[-1]
#     native=f'/gpfs/gpfs0/scratch/jws6pq/Notebook/PDB/3mer{protein_label}.pdb'
#     chi=f'/gpfs/gpfs0/scratch/jws6pq/Notebook/PDB/3merSARS2w{protein_label}S1.pdb'
#     try:
#         native_ssnf_index =ContactMap.correct_alignment_for_residue_position(aln,protein_label,index)
#         chimera_ssnf_index=ContactMap.get_sequence_from_pdb(chi,'C').find(soi)
#         native_contacts = ContactMap.get_individual_contacts(native, 'C', native_ssnf_index) + \
#                           ContactMap.intra_protein_contacts(native, 'C')[native_ssnf_index]
#         chi_contacts = ContactMap.get_individual_contacts(chi, 'C', chimera_ssnf_index) + \
#                        ContactMap.intra_protein_contacts(chi, 'C')[chimera_ssnf_index]
#         info=[protein_label,native_contacts,chi_contacts]
#     except:
#         info=[protein_label,'-','-']
#     print(info)
#     comparison_matrix[indexes]=info
# Analysis.convert_array_to_file(comparison_matrix,'/','/gpfs/gpfs0/scratch/jws6pq/Notebook/Fastas/Immersion/1792_contacts.csv')
# negative_index = [1266, 1265, 1264, 1030, 1267, 1263, 1270, 1262, 73, 1237, 1238, 1272, 1285, 1271, 1040, 1042, 749,
#                   856, 1287, 750, 1261, 1239, 751, 716, 1283, 1284, 501, 1286, 1396, 1297]

def convert_nucleotide_clusters_to_protein_alignment_index(cluster_file, columns_of_interest):
    nucleotide_protein_intersect = ()
    with open(cluster_file, 'r') as clusters:
        clusters = clusters.readlines()
    for cluster in clusters:
        cluster = cluster.split()
        corrected_cluster = [math.floor(int(index) / 3) for index in cluster]
        for index in corrected_cluster:
            if index in columns_of_interest:
                nucleotide_protein_intersect += (index,)
    return nucleotide_protein_intersect

def compare_aligned_contact(alignment_file, native_pdb, chimera_pdb, chain_id, label, alignment_index,sequence_of_interest,ref_label,
                            rank_difference_file=''):
    try:
        residue = ContactMap.get_residue_at_native_position(alignment_file, label, alignment_index)
        native_index = ContactMap.correct_alignment_for_residue_position(alignment_file, label, alignment_index)
    except:
        print('doesnt exist')
        return [label, alignment_index + 1, [], [], [], []]
    native_contacts = MDContacts.get_intra_residue_contact_list(native_pdb, chain_id,Path(native_pdb).stem+'_matrix')[native_index]
    native_seq = ContactMap.get_sequence_from_pdb(native_pdb, chain_id)
    if native_contacts:
        native_contacts_ids = [f'{native_seq[x]}{x + 1}' for x in native_contacts]
    else:
        native_contacts_ids = []
    native_contacts_ids += MDContacts.inter_residue_contact_list_md_analysis(native_pdb, chain_id, Path(native_pdb).stem+f'_{chain_id}_inter.json')[native_index]
    chimera_seq = ContactMap.get_sequence_from_pdb(chimera_pdb, chain_id)
    chi_index = ContactMap.correct_alignment_for_chimera_index(alignment_file, ref_label, alignment_index, label, sequence_of_interest)[0]
    residue_ids = f'{residue}{native_index + 1},{residue}{chi_index + 1}'
    chimera_contacts = MDContacts.get_intra_residue_contact_list(chimera_pdb, chain_id, Path(chimera_pdb).stem+'_matrix')[chi_index]
    if chimera_contacts:
        chimera_contacts_ids = [f'{chimera_seq[x]}{x + 1}' for x in chimera_contacts]
    else:
        chimera_contacts_ids = []
    chimera_contacts_ids += MDContacts.inter_residue_contact_list_md_analysis(chimera_pdb, chain_id, Path(chimera_pdb).stem+f'_{chain_id}_inter.json')[chi_index]
    # if rank_difference_file:
    #     rank_difference = rank_difference_table_to_dict(rank_difference_file)[label]
    #     return [label, alignment_index + 1, residue_ids, rank_difference, native_contacts_ids, chimera_contacts_ids]
    return [label, alignment_index + 1, residue_ids, native_contacts_ids, chimera_contacts_ids]

def compare_aligned_contact_all_proteins(alignment_file, native_pdb_format, chimera_pdb_format, chain_id,
                                         alignment_index, contact_comparison_file, sequence_of_interest,rank_difference_file=''):
    aln_dict=create_dictionary_from_alignment(alignment_file)
    comparison_matrix = zeros((len(aln_dict), 5 + int(bool(rank_difference_file))), dtype=object)

    for index, label in enumerate(aln_dict.keys()):
        native = native_pdb_format.format(label)
        chi = chimera_pdb_format.format(label)
        comparison = compare_aligned_contact(alignment_file, native, chi, chain_id, label, alignment_index,sequence_of_interest,
                                             '6VSB')
        comparison_matrix[index] = comparison
    Analysis.convert_array_to_file(comparison_matrix, '/', contact_comparison_file)


