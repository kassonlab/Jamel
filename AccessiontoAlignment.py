#! python
#
#  Code 2023 by Jamel Simpson

from os import system
import json
from pathlib import Path
from Bio import Entrez,Phylo,AlignIO,SeqIO
from random import choice
from Bio.Phylo.TreeConstruction import DistanceCalculator,DistanceTreeConstructor
from ChimeraGenerator import fasta_creation

def all_parents(tree):
    parents = dict()
    for node in tree.get_terminals():
        parents[tree.get_path(node)[-2]] = ()
    for node in tree.get_terminals():
        parents[tree.get_path(node)[-2]] += (node,)
    return parents


def parents_of_branch_terminals(tree,branch_clade):
    parents = dict()
    for node in tree.get_terminals():
        if branch_clade in tree.get_path(node):
            parents[tree.get_path(node)[-2]]=()
    for node in tree.get_terminals():
        if branch_clade in tree.get_path(node):
            parents[tree.get_path(node)[-2]] += (node,)
    return parents


def random_root_branch_children(newick_file,undesired_identifier=''):
    tree = Phylo.read(newick_file, 'newick')
    root_branches=[tree.get_path(node)[0] for node in tree.get_nonterminals() if len(tree.get_path(node))==1 ]
    terminal_parents={}
    if undesired_identifier:
        for branch in root_branches.copy():
            for terminal in branch.get_terminals():
                if terminal.name==undesired_identifier:
                    root_branches.remove(branch)
    for branch in root_branches:
        terminal_parents.update(parents_of_branch_terminals(tree,branch))
    diverse_selection=tuple(choice(children).name for parent,children in terminal_parents.items())
    return diverse_selection

def create_tree_from_aln(msa_file,new_tree_file='',new_tree_type='newick',new_ascii_file_representation=''):
    aln = AlignIO.read(msa_file, 'fasta')
    calculator = DistanceCalculator('identity')
    calculator.get_distance(aln)
    constructor = DistanceTreeConstructor(calculator, 'nj')
    tree = constructor.build_tree(aln)
    if new_tree_file:
        Phylo.write(tree,new_tree_file,new_tree_type)
    if new_ascii_file_representation:
        with open(new_ascii_file_representation, 'w') as handle:
            Phylo.draw_ascii(tree, file=handle)
    return tree


def accession_to_fasta(monomer_file_name, accession, email_for_Bio,subunits, multimer_name='NA'):
    """Takes an accession number and creates a fasta file with the sequence that corresponds with the accession given.
    A monomeric file is always created by default for alignment purposes even if a multimer file is requested"""
    Entrez.email = email_for_Bio
    # Pulling the sequence corresponding with accession numer specified
    handle = Entrez.efetch(db='protein', id=accession, retmode='text', rettype='fasta')
    # Turning the retrieved sequence into a single string with no breaks
    sequence = SeqIO.read(handle, "fasta").seq
    # Creating a monomer file by default for alignment purposes, if a multimer is requested it's made later
    fasta_creation(monomer_file_name, (sequence, 1,Path(monomer_file_name).stem))
    if subunits != 1:
        fasta_creation(multimer_name, (sequence, subunits,Path(multimer_name).stem))


def get_accession_sequence(accession, email_for_Bio):
    Entrez.email = email_for_Bio
    # Pulling the sequence corresponding with accession numer specified
    handle = Entrez.efetch(db='protein', id=accession, retmode='text', rettype='fasta')
    print(handle)
    return SeqIO.read(handle, "fasta").seq
# TODO make this fancy
def accession_to_fasta_nucleic(accession,  email_for_Bio,monomer_file_name=''):
    """Takes an accession number and creates a fasta file with the sequence that corresponds with the accession given.
    A monomeric file is always created by default for alignment purposes even if a multimer file is requested"""
    Entrez.email = email_for_Bio
    # Pulling the sequence corresponding with accession numer specified
    handle = Entrez.esearch(db='gene', term=accession,idtype='acc',rettype='uilist',retmode='json').readlines()[0]
    diction=json.loads(handle)

    try:
        gi = diction['esearchresult']['idlist'][0]
        handle = Entrez.efetch(db='gene', id=gi, retmode='text', rettype='seqid').readlines()
    except:
        print(monomer_file_name)
        return
    for x in handle:
        if 'Annotation' in x:
            nc_accession=x.split()[1]
            boundaries=x.split()[2].replace('(','').replace(')','').split('..')
            start=int(boundaries[0])-1
            end=int(boundaries[1])
            handle = Entrez.efetch(db='nuccore', id=nc_accession, retmode='text', rettype='fasta').readlines()
            sequence = ''.join(x for x in handle if x[0] != '>' if x != '').strip().replace('\n', '')[start:end]
            # Creating a monomer file by default for alignment purposes, if a multimer is requested it's made later
            if monomer_file_name:
                fasta_creation(monomer_file_name, (sequence, 1,Path(monomer_file_name).stem))
            return sequence


def create_dictionary_from_alignment(alignment_file):
    """Takes a fasta style alignment and makes a dictionary where the key is whatever signifier follows '>'
    and the value is the sequence with no spaces"""
    sequence_dictionary={}
    with open(alignment_file) as handle:
        for seq in SeqIO.parse(handle, "fasta"):

            sequence_dictionary[seq.id]=str(seq.seq)
    return sequence_dictionary



def multiple_sequence_fasta(list_of_sequences,new_fasta_file,list_of_fastas=''):

    if list_of_fastas:
        list_of_sequences=[]
        for fasta in list_of_fastas:
            list_of_sequences.append()

def multiple_sequence_alignment(list_of_fastas, fasta_for_alignment, new_alignment_file, reference_protein_fasta,muscle_command):
    """Creates a multiple sequence alignment using muscle and a concatenated fasta file with a reference fasta as the base,
     joined with all fastas specified in list_of_fastas."""
    # Creating a copy of the fasta file of a reference_protein_fasta to be added into
    # TODO use fasta creation formchimeragenerato
    system(f'cp  {reference_protein_fasta} {fasta_for_alignment}')
    # Concatenating the fasta files found in list_of_fasta
    for fasta in list_of_fastas:
        system(f"cat {fasta} >> {fasta_for_alignment}")
    # Using muscle to perform the alignment
    system(f'{muscle_command} -in {fasta_for_alignment} -out {new_alignment_file}')


def run_emboss_needle(new_emboss_file, sequence_one, sequence_two, needle_directory):
    """This runs EMBOSS on the command line."""
    system(f'{needle_directory}  -sprotein -gapopen 10 -gapextend 0.5 '
           f'-outfile {new_emboss_file} -asequence asis:{sequence_one} -bsequence asis:{sequence_two}')


def alignment_finder(alignment_file, sequence_of_interest, comparison_protein,
                     reference_protein, run_emboss='', new_emboss_file=''):
    """Takes a fasta style alignment and a sequence_of_interest from a reference_protein and returns the sequence of the
    comparison_protein that is outlined in the boundaries of the sequence_of_interest, as well as the python index boundaries
     for the found_alignment of the comparison_protein. reference_protein and comparison_protein must the names following
      '>' in the alignment file"""
    #Must be in fasta format
    with open(alignment_file, "r") as alignment:
        alignment = alignment.read().split('>')
        # Splitting up the sequences names and sequences into a dictionary
        sequence_dictionary = {sequence.split('\n')[0]: ''.join(sequence.split('\n')[1:]) for sequence in alignment if
                           len(sequence) != 0}
    # Recording the specified sequences as variables
    reference_sequence = sequence_dictionary[reference_protein]
    comparison_sequence = sequence_dictionary[comparison_protein]
    # Matching python indexing for the indexing from the alignment with some amount of '-' and indexing in the regular sequence
    reference_alignment_indexing = tuple(get_alignment_indexing(reference_sequence))
    # Creating a regular sequence without '-'
    no_gap_reference_sequence = ''.join(x for ind, x in enumerate(reference_sequence) if x.isalpha())
    # Boundaries are given in python index
    alignment_reference_start = reference_alignment_indexing[no_gap_reference_sequence.find(sequence_of_interest)]
    alignment_reference_end = reference_alignment_indexing[no_gap_reference_sequence.find(sequence_of_interest) + len(sequence_of_interest)-1]
    # Pulling the section of the comparison_sequence that overlaps with the sequence_of_interest
    found_alignment = comparison_sequence[alignment_reference_start:(alignment_reference_end+1)].replace('-', '')
    no_gap_reference_start = no_gap_reference_sequence.find(sequence_of_interest)
    no_gap_reference_end = no_gap_reference_start + len(sequence_of_interest)
    # Recording the indexes of the found_alignment
    # Additionally the splice_start is the first residue that is spliced,
    # and splice_end is the first residue that's not spliced
    splice_start = comparison_sequence.replace('-', '').find(found_alignment)
    splice_end = splice_start + len(found_alignment)
    if run_emboss != '':
        run_emboss_needle(new_emboss_file, found_alignment, sequence_of_interest, run_emboss)
    return found_alignment, (splice_start, splice_end), (no_gap_reference_start, no_gap_reference_end)


def truncate_alignmnent(alignment_file,sequence_of_interest,new_aln_file):
    """Might be dilapidated"""
    with open(alignment_file, 'r') as alignment:
        alignment = alignment.read().split('>')
        sequence_dictionary = {sequence.split('\n')[0]: ''.join(sequence.split('\n')[1:]) for sequence in alignment
                                        if
                                        len(sequence) != 0}
        reference_sequence = sequence_dictionary['6vsb_adjacent']
        # Matching python indexing for the indexing from the alignment with some amount of '-' and indexing in the regular sequence
        reference_alignment_indexing = tuple((ind for ind, x in enumerate(reference_sequence) if x.isalpha()))
        no_gap_reference_sequence = ''.join(x for ind, x in enumerate(reference_sequence) if x.isalpha())
        alignment_reference_start = reference_alignment_indexing[no_gap_reference_sequence.find(sequence_of_interest)]
        alignment_reference_end = reference_alignment_indexing[
                                      no_gap_reference_sequence.find(sequence_of_interest) + len(sequence_of_interest) - 1] + 1
        list_of_tuples=[]
        for stem,sequence in sequence_dictionary.items():
            list_of_tuples.append((sequence[alignment_reference_start:alignment_reference_end],1,stem))
            fasta_creation(new_aln_file,list_of_tuples)


def check_mutation_cutoff(cutoff,seq_1,seq_2):
    difference_count=0
    for res1,res2 in zip(seq_1,seq_2):
        difference_count+=int(res1!=res2)
        if difference_count>=cutoff:
            return True
    return False

def exclude_related_sequences(cutoff,alignment_file,starting_sequence,included_file='',excluded_file=''):
    with open(alignment_file, "r") as alignment:
        alignment = alignment.read().split('>')
        sequence_dictionary = {sequence.split('\n')[0]: ''.join(sequence.split('\n')[1:]) for sequence in alignment
                           if
                           len(sequence) != 0}
    included = {starting_sequence:sequence_dictionary[starting_sequence]}
    excluded={}
    for strain, sequence in sequence_dictionary.items():
        sufficiently_distant=0
        for accepted_sequence in included.copy().values():
            if not check_mutation_cutoff(cutoff,sequence,accepted_sequence):
                break
            sufficiently_distant+=1
        if sufficiently_distant==len(included.keys()):
            included[strain] = sequence
        else:
            excluded[strain]=sequence
    print('excluded:',len(excluded))
    print('included:', len(included))
    if included_file:
        with open(included_file, 'w') as outfile:
            for strain,sequence in included.items():
                outfile.write(f'>{strain}\n{sequence}\n')
    if excluded_file:
        with open(excluded_file, 'w') as outfile:
            for strain,sequence in excluded.items():
                outfile.write(f'>{strain}\n{sequence}\n')
def create_list_of_fasta_files(list_of_fastas,file_name):
    with open(file_name, 'w') as fasta_list_file:
        fasta_list_file.write("\n".join(list_of_fastas))

def get_alignment_indexing(alignment_seq):
    """Creating a list of alignment indexes for non '-' characters"""
    return [ind for ind, x in enumerate(alignment_seq) if x.isalpha()]
def no_gap_sequence_from_alignment(alignment_seq):
    """Removes gaps from an alignment sequence"""
    return alignment_seq.replace('-','')

def clustalw_to_fasta(clustal_aln_file,new_fasta_aln_file):
    """Converts clustal w aliignment file into a fasta alignment"""
    clw = SeqIO.parse(clustal_aln_file, "clustal")
    return SeqIO.write(clw, new_fasta_aln_file, "fasta")

# ODO needs to be a monomer
def extract_seq_from_fasta(fasta_file):
    with open(fasta_file, "r") as handle:
        return SeqIO.read(handle, "fasta").seq

print(extract_seq_from_fasta('6vsb_S1.fasta'))
# def align_pdb_sequences(ref_pdb_chain_tuple,comparison_pdb_chain_tuple,new_alignment,muscle_command):
#     get_sequence_from_pdb()


