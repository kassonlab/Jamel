#! python
#
#  Code 2023 by Jamel Simpson
import re
import subprocess
from os import system
import json
from pathlib import Path
from Bio import Entrez, Phylo, AlignIO, SeqIO
from Bio.Align import PairwiseAligner
from random import choice
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def all_parents(tree):
    parents = dict()
    for node in tree.get_terminals():
        parents[tree.get_path(node)[-2]] = ()
    for node in tree.get_terminals():
        parents[tree.get_path(node)[-2]] += (node,)
    return parents


def parents_of_branch_terminals(tree, branch_clade):
    parents = dict()
    for node in tree.get_terminals():
        if branch_clade in tree.get_path(node):
            parents[tree.get_path(node)[-2]] = ()
    for node in tree.get_terminals():
        if branch_clade in tree.get_path(node):
            parents[tree.get_path(node)[-2]] += (node,)
    return parents


def random_root_branch_children(newick_file, undesired_identifier=''):
    tree = Phylo.read(newick_file, 'newick')
    root_branches = [tree.get_path(node)[0] for node in tree.get_nonterminals() if len(tree.get_path(node)) == 1]
    terminal_parents = {}
    if undesired_identifier:
        for branch in root_branches.copy():
            for terminal in branch.get_terminals():
                if terminal.name == undesired_identifier:
                    root_branches.remove(branch)
    for branch in root_branches:
        terminal_parents.update(parents_of_branch_terminals(tree, branch))
    diverse_selection = tuple(choice(children).name for parent, children in terminal_parents.items())
    return diverse_selection


def create_tree_from_aln(msa_file, new_tree_file='', new_tree_type='newick', new_ascii_file_representation=''):
    aln = AlignIO.read(msa_file, 'fasta')
    calculator = DistanceCalculator('identity')
    calculator.get_distance(aln)
    constructor = DistanceTreeConstructor(calculator, 'nj')
    tree = constructor.build_tree(aln)
    if new_tree_file:
        Phylo.write(tree, new_tree_file, new_tree_type)
    if new_ascii_file_representation:
        with open(new_ascii_file_representation, 'w') as handle:
            Phylo.draw_ascii(tree, file=handle)
    return tree


def translate_dna_to_protein(dna_seq):
    return Seq(dna_seq).translate()


def fasta_creation(file_name, sequences: list[SeqRecord]):
    """Creates a fasta file with the given file_name, and replicates the sequence within it the specified number of times
    to create a homo multimer if subunits is greater than 1."""
    with open(file_name, 'w') as outfile:
        SeqIO.write(sequences, outfile, "fasta")


def accession_to_fasta(accession, email_for_bio, subunits, fasta_id):
    """Takes an accession number and creates a fasta file with the sequence that corresponds with the accession given.
    A monomeric file is always created by default for alignment purposes even if a multimer file is requested"""
    Entrez.email = email_for_bio
    # Pulling the sequence corresponding with accession numer specified
    handle = Entrez.efetch(db='protein', id=accession, retmode='text', rettype='fasta')
    # Turning the retrieved sequence into a single string with no breaks
    sequence = SeqIO.read(handle, "fasta").seq
    # Creating a monomer file by default for alignment purposes, if a multimer is requested it's made later
    fasta_creation(fasta_id + '.fa', [SeqRecord(Seq(sequence), id=fasta_id, description='') for x in range(subunits)])


def create_seq_records(label, seq, description='', subunit_count=1):
    return [SeqRecord(Seq(seq), id=label, description=description) for _ in range(subunit_count)]


def get_accession_sequence(accession, email_for_bio):
    Entrez.email = email_for_bio
    # Pulling the sequence corresponding with accession numer specified
    handle = Entrez.efetch(db='protein', id=accession, retmode='text', rettype='fasta')
    return SeqIO.read(handle, "fasta").seq


def accession_to_fasta_nucleic(accession, email_for_Bio, monomer_file_name=''):
    """Takes an accession number and creates a fasta file with the sequence that corresponds with the accession given.
    A monomeric file is always created by default for alignment purposes even if a multimer file is requested"""
    Entrez.email = email_for_Bio
    # Pulling the sequence corresponding with accession numer specified
    handle = Entrez.esearch(db='gene', term=accession, idtype='acc', rettype='uilist', retmode='json').readlines()[0]
    diction = json.loads(handle)

    try:
        gi = diction['esearchresult']['idlist'][0]
        handle = Entrez.efetch(db='gene', id=gi, retmode='text', rettype='seqid').readlines()
    except:
        print(monomer_file_name)
        return None
    for x in handle:
        if 'Annotation' in x:
            nc_accession = x.split()[1]
            boundaries = x.split()[2].replace('(', '').replace(')', '').split('..')
            start = int(boundaries[0]) - 1
            end = int(boundaries[1])
            handle = Entrez.efetch(db='nuccore', id=nc_accession, retmode='text', rettype='fasta').readlines()
            sequence = ''.join(x for x in handle if x[0] != '>' if x != '').strip().replace('\n', '')[start:end]
            # Creating a monomer file by default for alignment purposes, if a multimer is requested it's made later
            if monomer_file_name:
                fasta_creation(monomer_file_name,
                               create_seq_records(Path(monomer_file_name).stem, sequence, subunit_count=1))
            return sequence
        return None
    return None


def create_dictionary_from_alignment(alignment_file):
    """Takes a fasta style alignment and makes a dictionary where the key is whatever signifier follows '>'
    and the value is the sequence with no spaces"""
    with open(alignment_file) as handle:
        return {seq.id: str(seq.seq) for seq in SeqIO.parse(handle, "fasta")}


def dictionary_to_fasta(seq_dict: dict[str:str], new_fasta_file):
    seq_records = [SeqRecord(Seq(sequence), id=seq_id, description="") for seq_id, sequence in seq_dict.items()]
    with open(new_fasta_file, "w") as output_handle:
        SeqIO.write(seq_records, output_handle, "fasta")


def multiple_sequence_alignment(sequences: list[SeqRecord], fasta_for_alignment, new_alignment_file, muscle_command):
    """Creates a multiple sequence alignment using muscle and a concatenated fasta file with a reference fasta as the base,
     joined with all fastas specified in list_of_fastas."""
    # Creating a copy of the fasta file of a reference_protein_fasta to be added into
    fasta_creation(fasta_for_alignment, sequences)
    # Using muscle to perform the alignment
    system(f'{muscle_command} -align {fasta_for_alignment} -output {new_alignment_file}')


def run_emboss_needle(new_emboss_file, sequence_one: str, sequence_two: str, emboss_command='needle'):
    """This runs EMBOSS on the command line."""
    subprocess.run(f'{emboss_command}  -sprotein -gapopen 10 -gapextend 0.5 '
                   f'-outfile {new_emboss_file} -asequence asis:{sequence_one} -bsequence asis:{sequence_two}',
                   shell=True)


def calculate_sequence_identity(aln_seq_1='', aln_seq_2='', seq_1='', seq_2=''):
    if seq_2 and seq_1:
        aligner = PairwiseAligner()
        aligner.mode = 'global'
        return calculate_sequence_identity(*aligner.align(seq_1, seq_2)[0])
    elif aln_seq_1 and aln_seq_2:
        if len(aln_seq_1) != len(aln_seq_2): raise ValueError(
            f'Alignment sequences are not equal length. aln_seq_1:{len(aln_seq_1)} aln_seq_2:{len(aln_seq_2)}')
        identical_count = [resb == resa for resa, resb in zip(aln_seq_1, aln_seq_2) if resa != '-' or resb != '-']
        return identical_count.count(True) / len(identical_count)
    else:
        raise ValueError('Did not provide sufficient pairs of either aligned or unligned sequences.')


def get_alignment_indexing(alignment_seq):
    """Creating a list of indexes that represent the positions of non-gap amino acids in an alignment"""
    return [ind for ind, x in enumerate(alignment_seq) if x != '-']


def get_alignment_indexing_w_dashes(alignment_seq):
    return [ind if x != '-' else '-' for ind, x in enumerate(alignment_seq)]


def no_gap_sequence_from_alignment(alignment_seq):
    """Removes gaps from an alignment sequence"""
    return re.sub(r'[^a-zA-Z]', '', alignment_seq)


def find_overlapping_matches(string, substring):
    return [m.span()[0] for m in re.finditer(rf'(?={substring})', string)]


def alignment_finder(sequence_of_interest, partner_label, base_label, aln_file):
    """Takes a fasta style alignment and a sequence_of_interest from a base_protein and returns the sequence of the
    comparison_protein that aligns with the sequence_of_interest in the base sequence, as well as the python index boundaries
     for the found_seq of the comparison_protein. base_protein and comparison_protein must be the names following
      '>' in the alignment file"""
    # Matching python indexing for the indexing from the alignment with some amount of '-' and indexing in the regular sequence
    aln = create_dictionary_from_alignment(aln_file)
    base_aln = aln[base_label]
    partner_aln = aln[partner_label]
    base_aln_indexing = get_alignment_indexing(base_aln)
    # Creating a regular sequence without '-'
    base_seq = no_gap_sequence_from_alignment(base_aln)
    interest_len = len(sequence_of_interest)
    matches = find_overlapping_matches(base_seq, sequence_of_interest)
    if len(matches) == 0:
        raise ValueError('Could not find sequence of interest in base sequence')
    if len(matches) > 1:
        raise ValueError('Sequence of interest is non-specific. Found multiple occurences')
    else:
        base_start = matches[0]
        base_end = base_start + interest_len - 1
        # Boundaries are given in python index
        alignment_base_start = base_aln_indexing[base_start]
        # Because indexing is exclusive the final index is put up one
        alignment_base_end = base_aln_indexing[base_end] + 1
        # Pulling the section of the comparison_sequence that overlaps with the sequence_of_interest
        found_seq = no_gap_sequence_from_alignment(partner_aln[alignment_base_start:alignment_base_end])
        if found_seq:
            # I'm finding the first occurence of an amino acid in the alignment that overlaps with the sequence of interest
            # in order to calculate the indices in the regular sequence
            partner_indexing = get_alignment_indexing(partner_aln)
            min_abs = sorted(partner_indexing, key=lambda x: abs(alignment_base_start - x))
            lowest_2 = list(alignment_base_start - x for x in min_abs[0:2])
            if abs(lowest_2[0]) == abs(lowest_2[1]):
                closest_aln_pos = min(min_abs[0:2], key=lambda x: alignment_base_start - x)
            else:
                closest_aln_pos = min_abs[0]
            partner_start = partner_indexing.index(closest_aln_pos)
            partner_end = partner_start + len(found_seq) - 1
        else:
            # if there are no amino acids in the complement (meaning found_seq is empty) then the inheritance dict will show accordingly
            partner_start = partner_end = 0
        chimera_aln_seq = base_aln.replace(base_aln[alignment_base_start:alignment_base_end],
                                           partner_aln[alignment_base_start:alignment_base_end])
        # this inheritance dictionary shows the splice boundaries from a given parent and gives the corresponding boundary in the chimera where that parent sequence is.
        # data structure is {parent:{parent_splice:chimera_splice}}
        inheritance = make_inheritance_dict(base_label, base_seq, (base_start, base_end), partner_label,
                                            (partner_start, partner_end),
                                            no_gap_sequence_from_alignment(chimera_aln_seq), len(found_seq))
        splice_boundaries={'base_start':base_start,'base_end':base_end,'partner_start':partner_start,'partner_end':partner_end}
        return chimera_aln_seq, inheritance,splice_boundaries


def make_inheritance_dict(base_label: str, base_seq: str, base_boundaries: tuple, partner_label: str,
                          partner_boundaries: tuple, chimera_seq: str, insertion_length: int):
    base_start, base_end = base_boundaries
    partner_start, partner_end = partner_boundaries
    return {base_label: {(0, base_start): (0, base_start), (base_end, len(base_seq)): (
        base_start + insertion_length - 1, len(chimera_seq))},
            partner_label: {(partner_start, partner_end): (base_start, base_start + insertion_length - 1)}}

def convert_fasta_aln_to_clustal(fasta_file,new_clustal):
    records = SeqIO.parse(fasta_file, "fasta")
    count = SeqIO.write(records, new_clustal, "clustal")

def map_plddt_to_aln(aln_seq, plddt):
    aln_index = get_alignment_indexing(aln_seq)
    aln_plddt = []
    for aln_pos, res in enumerate(aln_seq):
        score = plddt[aln_index.index(aln_pos)] if aln_pos in aln_index else res
        aln_plddt.append(score)
    return aln_plddt


def contiguous_inheritance_dict(chi_label: str, parent_labels: tuple, aln_file: str):
    """Mostly designed for the SCHEMA dataset, It finds all overlapping amino acids between parents sequences of a chimera,
    as opposed to the inheritance dict from block_swap_inheritance that is used for single continuous block splices"""
    residue_inheritance = {parent: set() for parent in parent_labels}
    seq_dict = create_dictionary_from_alignment(aln_file)
    for parent in parent_labels:
        chi_aln = seq_dict[chi_label]
        residue_inheritance[parent] = {aln_pos for aln_pos, (chi_res, par_res) in enumerate(
            zip(chi_aln, seq_dict[parent])) if chi_res == par_res if chi_res != '-'}
    return residue_inheritance


def check_mutation_cutoff(cutoff, seq_1, seq_2):
    difference_count = 0
    for res1, res2 in zip(seq_1, seq_2):
        difference_count += int(res1 != res2)
        if difference_count >= cutoff:
            return True
    return False


def exclude_related_sequences(cutoff, alignment_file, starting_sequence, included_file='', excluded_file=''):
    """Filters alignment sequences based on provided sequence identity cutoff"""
    seq_dict = create_dictionary_from_alignment(alignment_file)
    included = {starting_sequence: seq_dict[starting_sequence]}
    excluded = {}
    for strain, sequence in seq_dict.items():
        sufficiently_distant = 0
        for accepted_sequence in included.copy().values():
            if not check_mutation_cutoff(cutoff, sequence, accepted_sequence):
                break
            sufficiently_distant += 1
        if sufficiently_distant == len(included.keys()):
            included[strain] = sequence
        else:
            excluded[strain] = sequence
    print('excluded:', len(excluded))
    print('included:', len(included))
    if included_file:
        dictionary_to_fasta(included, included_file)
    if excluded_file:
        dictionary_to_fasta(excluded, excluded_file)


def create_list_of_fasta_files(list_of_fastas, file_name):
    with open(file_name, 'w') as fasta_list_file:
        fasta_list_file.write("\n".join(list_of_fastas))


def clustalw_to_fasta(clustal_aln_file, new_fasta_aln_file):
    """Converts clustal w aliignment file into a fasta alignment"""
    clw = SeqIO.parse(clustal_aln_file, "clustal")
    return SeqIO.write(clw, new_fasta_aln_file, "fasta")


def extract_seq_from_fasta(fasta_file):
    with open(fasta_file, "r") as handle:
        return str(SeqIO.read(handle, "fasta").seq)


def straighten_alignment(aln_file, new_aln_file):
    with open(aln_file, "r") as aln, open(new_aln_file, "w") as out:
        records = SeqIO.parse(aln, "fasta")
        SeqIO.write(records, out, 'fasta-2line')
