import math
from Bio import Seq
import Analysis
from AccessiontoAlignment import get_alignment_indexing, no_gap_sequence_from_alignment, clustalw_to_fasta, \
    alignment_finder,create_dictionary_from_alignment
from numpy import zeros, savetxt
from scipy.stats import rankdata
import MDContacts
import  ContactMap

'''This is very specific code for dilapidated and wrong sectioning of S1 and conversion of individual alignments to 
map on a MSA'''

vsb_S1 = 'QCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPGSASS'


def dilapidated_translate_old_to_new_alignment(old_alignment_file, reference_label, comparison_label, alignment_file,
                                               sequence_of_interest):
    """This is a dilapidated piece of code to find the boundaries of overlapping sequences in a pairwise fasta
    alignment, but it has a purposeful mistake from previous iterations of alignment_finder that will help replicate previous chimeric
    sequences until new data with more accurate chimeras are created. Takes a fasta style alignment and a
    sequence_of_interest from a reference_protein and returns the sequence of the comparison_protein that is outlined
    in the boundaries of the sequence_of_interest, as well as the python index boundaries for the found_alignment of
    the comparison_protein. reference_protein and comparison_protein must the names following '>' in the alignment
    file"""
    old_sequence_dict = create_dictionary_from_alignment(old_alignment_file)
    old_reference_alignment = old_sequence_dict[reference_label]
    old_comparison_alignment = old_sequence_dict[comparison_label]
    reference_alignment_indexing = get_alignment_indexing(old_reference_alignment)
    no_gap_reference_sequence = no_gap_sequence_from_alignment(old_reference_alignment)
    # Index where the sequence of interest starts
    alignment_reference_start = reference_alignment_indexing[no_gap_reference_sequence.find(sequence_of_interest)]
    alignment_reference_end = reference_alignment_indexing[
        no_gap_reference_sequence.find(sequence_of_interest) + len(sequence_of_interest)]
    # This is the mistake, where the erasure of gaps happens prematurely before the correct sequence occurs. However,
    # this generates what is the most recent chimeric slice for the data as of 07/2023. It should read
    # 'old_comparison_sequence[alignment_reference_start:alignment_reference_end].replace('-', '')'
    old_sequence_chunk = old_comparison_alignment.replace('-', '')[alignment_reference_start:alignment_reference_end]
    sequence_dict = create_dictionary_from_alignment(alignment_file)
    comparison_alignment = sequence_dict[comparison_label]
    comparison_alignment_indexing = get_alignment_indexing(comparison_alignment)
    comparison_sequence = no_gap_sequence_from_alignment(comparison_alignment)
    new_alignment_comparison_start = comparison_alignment_indexing[comparison_sequence.find(old_sequence_chunk)]
    new_alignment_comparison_end = comparison_alignment_indexing[
        comparison_sequence.find(old_sequence_chunk) + len(old_sequence_chunk)]
    return old_sequence_chunk, new_alignment_comparison_start, new_alignment_comparison_end


def alignment_to_confidence(aln_file, comparison_label, ref_label, pdb_file, soi):
    '''This only works with homologous proteins for now'''
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
    confidence_scores = listed_alignment_characters
    return confidence_scores


# turn this into smaller so you can put soi for each indivdual
def confidence_rank_matrix(alignment_file, pdb_naming, reference_label, rank_matrix_file, sequence_of_interest,
                           raw_matrix_file=''):
    aln_dict = create_dictionary_from_alignment(alignment_file)
    aln_seq_length = len(aln_dict[reference_label])
    aln_length = len(aln_dict)
    plddt_matrix = zeros((aln_length + 1, aln_seq_length + 1), dtype=object)
    for index, label in enumerate(aln_dict.keys()):
        plddt_matrix[index, 0] = label
        plddt_matrix[index, 1:] = alignment_to_confidence(alignment_file, label, reference_label,
                                                          pdb_naming.replace('*', label), sequence_of_interest)
    if raw_matrix_file:
        savetxt(raw_matrix_file, plddt_matrix, fmt='%s')
    for index, column in enumerate(plddt_matrix[0, 1:]):
        column = plddt_matrix[:, 1 + index]
        rank_column = rankdata([x for x in column if x != '-'], 'ordinal')
        column_indexing = [ind for ind, x in enumerate(column) if x != '-']
        for rank_index, correct_index in enumerate(column_indexing):
            column[correct_index] = rank_column[rank_index]
        plddt_matrix[:, 1 + index] = column
    savetxt(rank_matrix_file, plddt_matrix, fmt='%s')


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


def contact_contingency(alignment_file, native_pdb, chimera_pdb, chain_id, label, alignment_index, rank_difference_file,
                        contact_to_find_aln_index,ref_label,sequence_of_interest):
    try:
        residue = ContactMap.get_residue_at_native_position(alignment_file, label, alignment_index)
        native_index = ContactMap.correct_alignment_for_residue_position(alignment_file, label, alignment_index)
    except:
        print('doesnt exist')
        return [label, '-', '-', '-', '-', '-', '-']
    native_contacts = MDContacts.intra_residue_dist_matrix_md_analysis(native_pdb, chain_id)[native_index]
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

    chi_index = ContactMap.correct_alignment_for_chimera_index(alignment_file, ref_label, alignment_index, label,sequence_of_interest)
    residue_ids = f'{residue}{native_index + 1},{residue}{chi_index + 1}'
    chimera_contacts = ContactMap.get_individual_intra_contacts(chimera_pdb, chain_id, chi_index)
    chimera_contacts += ContactMap.get_inter_protein_contacts(chimera_pdb, chain_id)[chi_index]
    chimera_contacts = [int(str(contacts).split(':')[-1]) for contacts in chimera_contacts]
    chi_index_to_find = ContactMap.correct_alignment_for_chimera_index(alignment_file, ref_label, contact_to_find_aln_index, label,sequence_of_interest)
    chi_yes = chi_index_to_find in chimera_contacts
    chi_no = chi_index_to_find not in chimera_contacts
    rank_difference = confidence_rank_matrix(rank_difference_file)[label]
    return [label, residue_ids, rank_difference, int(native_yes), int(native_no), int(chi_yes), int(chi_no),
            residue_to_find]


# with open("/gpfs/gpfs0/scratch/jws6pq/Notebook/Overall/List_of_coronaviruses", 'r') as loc:
#     loc = loc.readlines()
# aln='/gpfs/gpfs0/scratch/jws6pq/Notebook/PDB/CoronavirusMSA.aln'
# comparison_matrix=zeros((len(loc),8),dtype=object)
# index=1267
# rank_change = f'/gpfs/gpfs0/scratch/jws6pq/Notebook/PDB/Rank_change_{index+1}.tsv'
# for indexes,label in enumerate(loc):
#     protein_label=label.split()[-1]
#     native=f'/gpfs/gpfs0/scratch/jws6pq/Notebook/PDB/3mer{protein_label}.pdb'
#     chi=f'/gpfs/gpfs0/scratch/jws6pq/Notebook/PDB/3merSARS2w{protein_label}S1.pdb'
#     comparison=contact_contingency(aln, native, chi, 'B', protein_label,index, rank_change,1792,'SSNF')
#     print(comparison)
#     comparison_matrix[indexes]=comparison
# savetxt(f'/gpfs/gpfs0/scratch/jws6pq/Notebook/Immersion/CarbonB_1268_contacts.csv', comparison_matrix, fmt='%s/%s/%s/%s/%s/%s/%s/%s')

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



