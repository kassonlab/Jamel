import Analysis
from Bio import SeqIO

from numpy import zeros,savetxt
from scipy.stats import rankdata
'''This is very specific code for dilapidated and wrong sectioning of S1 and conversion of individual alignments to map on a MSA'''
def alignment_to_confidence(old_alignment,reference_label, comparison_label, alignment_file, sequence_of_interest, pdb):
    with open(old_alignment, "r") as alignment:
        alignment = alignment.read().split('>')
        # Splitting up the sequences names and sequences into a dictionary
        sequence_dictionary = {sequence.split('\n')[0]: ''.join(sequence.split('\n')[1:]) for sequence in alignment if
                               len(sequence) != 0}
    # Recording the specified sequences as variables
    reference_sequence = sequence_dictionary[reference_label]
    comparison_sequence = sequence_dictionary[comparison_label]
    # Matching python indexing for the indexing from the alignment with some amount of '-' and indexing in the regular sequence
    reference_alignment_indexing = tuple((ind for ind, x in enumerate(reference_sequence) if x.isalpha()))
    # TODO look at old alignment finders around late 2022 in google chrome
    # Creating a regular sequence without '-'
    no_gap_reference_sequence = ''.join(x for ind, x in enumerate(reference_sequence) if x.isalpha())
    # Boundaries are given in python index
    alignment_reference_start = reference_alignment_indexing[no_gap_reference_sequence.find(sequence_of_interest)]
    alignment_reference_end = reference_alignment_indexing[
        no_gap_reference_sequence.find(sequence_of_interest) + len(sequence_of_interest)]
    old_sequence_chunk=comparison_sequence.replace('-','')[alignment_reference_start:alignment_reference_end]
    with open(alignment_file, "r") as alignment:
        alignment = alignment.read().split('>')
        # Splitting up the sequences names and sequences into a dictionary
        sequence_dictionary = {sequence.split('\n')[0]: ''.join(sequence.split('\n')[1:]) for sequence in alignment if
                               len(sequence) != 0}
    # Recording the specified sequences as variables
    comparison_sequence = sequence_dictionary[comparison_label]
    # Matching python indexing for the indexing from the alignment with some amount of '-' and indexing in the regular sequence
    comparison_alignment_indexing = tuple((ind for ind, x in enumerate(comparison_sequence) if x.isalpha()))
    # Creating a regular sequence without '-'
    no_gap_comparison_sequence = ''.join(x for ind, x in enumerate(comparison_sequence) if x.isalpha())
    # Boundaries are given in python index
    alignment_comparison_start = comparison_alignment_indexing[no_gap_comparison_sequence.find(old_sequence_chunk)]
    alignment_comparison_end = comparison_alignment_indexing[
        no_gap_comparison_sequence.find(old_sequence_chunk) + len(old_sequence_chunk)]
    # Pulling the section of the comparison_sequence that overlaps with the sequence_of_interest
    total_alignment_list=list(comparison_sequence)
    found_alignment_indexing = comparison_alignment_indexing[no_gap_comparison_sequence.find(old_sequence_chunk):no_gap_comparison_sequence.find(old_sequence_chunk) + len(old_sequence_chunk)]
    no_gap_comparison = comparison_sequence[alignment_comparison_start:alignment_comparison_end].replace('-', '')
    plddt = Analysis.get_plddt_dict_from_pdb(pdb)
    for sequence, plddt in plddt.items():
        start = sequence.find(no_gap_comparison)
        end = start + len(no_gap_comparison)
        plddt_chunk = plddt[start:end]
        for index, position in enumerate(found_alignment_indexing):
            total_alignment_list[position] = plddt_chunk[index]
    return total_alignment_list

def confidence_rank_matrix(alignment_file, list_of_label_pdb_tuples, reference_label, rank_matrix_file,
                           raw_matrix_file=''):
    with open(alignment_file, "r") as alignment:
        alignment = alignment.read().split('>')
        # Splitting up the sequences names and sequences into a dictionary
        alignment_length = len({sequence.split('\n')[0]: ''.join(sequence.split('\n')[1:]) for sequence in alignment if
                               len(sequence) != 0}[reference_label])
    plddt_matrix = zeros((len(list_of_label_pdb_tuples) + 1, alignment_length + 1), dtype=object)
    for index, (label, pdb) in enumerate(list_of_label_pdb_tuples):
        plddt_matrix[index, 0] = label
        print(label)
        plddt_matrix[index, 1:] = alignment_to_confidence(f"/gpfs/gpfs0/scratch/jws6pq/Notebook/Alignment/{label}onSARS2.fasta",reference_label, label, alignment_file, S1.replace('-', ''),
                                                          pdb)
    print('SARS')
    plddt_matrix[-1, 0] = 'SARS2'
    plddt_matrix[-1, 1:] = Analysis.alignment_to_confidence(reference_label, reference_label, alignment_file,
                                                   S1.replace('-', ''),
                                                   f'/gpfs/gpfs0/scratch/jws6pq/Notebook/PDB/3merSARS2.pdb')
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
S1 = '--------------------AYTNS-------------------------------------------------------------------------------------------------FTRGV-------YYPDK--------------------------VFR--------------------------------------------S-------SVLHSTQD----LFLPFFS-----------------------------------------------------NVTW-----FHNP----------------------VL--PFNDGVYFASTEKSNII--------------------------------RGWIFG---TTLD----------------SKTQSLLIVN-------------------------------NA--TN---------------------VVIKVC--------------------------EFQFC--NDP-------FLSE--FRVYSS------------------------------ANNCTFEYVSQPFLKNLRE--------------FVFKNIDG-YFKIYSK-------------------HTPPQGFSALEPLVDL-PIGINIT---RFQTLLAA------------------------YYVGYL------QPRTFLLKYNENGTITDAVDCALD-PLSETKCTLKSFTVEKGIYQTSNFRV--QPTESIV--RFPNITNLC--PFGEVFNAT--RFAS---VYAWNRKRIS--NC-----VADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSF--VIRGDEVRQIA--PGQTGKIADYNYKLPDDFTGCVIAWNSNNL-----------DSYNYL---------YRN-------------------LKP------------FER-----DIS----T-EIY--------------------------------------------------------------NC---YFP----------LQ---------------------SYGFQ------PT----VGYQPYR------------VVVLSFELL------------------------------------------------------HAPATVCGPK-----KS------------------TNLVKNKCVNFNFNGLTGTGVLTE-SN----KKFL-PFQQFGRDIAD--TTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVT-------------------------------GSNVFQTR----AGCLIGAEHVNNSY------E--CDIPI----GAGICA'
with open("/gpfs/gpfs0/scratch/jws6pq/Notebook/Overall/List_of_coronaviruses", 'r') as loc:
    loc = loc.readlines()
    label_pdb_list = []
    chi_label_pdb_list = []
    for line in loc:
        file_stem = line.split()[-1]
        label_pdb_list.append((file_stem, f'/gpfs/gpfs0/scratch/jws6pq/Notebook/PDB/3mer{file_stem}.pdb'))
        chi_label_pdb_list.append((file_stem, f'/gpfs/gpfs0/scratch/jws6pq/Notebook/PDB/3merSARS2w{file_stem}S1.pdb'))

        # records = SeqIO.parse(f"/gpfs/gpfs0/scratch/jws6pq/Notebook/Alignment/{file_stem}onSARS2.aln", "clustal")
        # count = SeqIO.write(records, f"/gpfs/gpfs0/scratch/jws6pq/Notebook/Alignment/{file_stem}onSARS2.fasta", "fasta")

confidence_rank_matrix("/gpfs/gpfs0/scratch/jws6pq/Notebook/Fastas/SARS_S1_chimeras/CoronavirusMSA.aln",
                       chi_label_pdb_list, '6vsb_B',
                       '/gpfs/gpfs0/scratch/jws6pq/Notebook/Fastas/Immersion/rank_chi_plddt_matrix.csv',
                       '/gpfs/gpfs0/scratch/jws6pq/Notebook/Fastas/Immersion/raw_chi_plddt_matrix.csv')
confidence_rank_matrix("/gpfs/gpfs0/scratch/jws6pq/Notebook/Fastas/SARS_S1_chimeras/CoronavirusMSA.aln",
                       label_pdb_list, '6vsb_B',
                       '/gpfs/gpfs0/scratch/jws6pq/Notebook/Fastas/Immersion/rank_native_plddt_matrix.csv',
                       '/gpfs/gpfs0/scratch/jws6pq/Notebook/Fastas/Immersion/raw_native_plddt_matrix.csv')