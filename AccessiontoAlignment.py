def accession_to_fasta(monomer_file_name, accession, email_for_Bio,subunits=3, multimer_name='NA'):
    """Takes an accession number and creates a fasta file with the sequence that corresponds with the accession given.
    A monomeric file is always created by default for alignment purposes even if a multimer file is requested"""
    from Bio import Entrez
    from ChimeraGenerator import fasta_creation
    Entrez.email = email_for_Bio
    # Pulling the sequence corresponding with accession numer specified
    handle = Entrez.efetch(db='protein', id=accession, retmode='text', rettype='fasta').readlines()
    # Turning the retrieved sequence into a single string with no breaks
    sequence=''.join([x for x in handle if  x[0]!='>'  if x!='']).strip().replace('\n','')
    # Creating a monomer file by default for alignment purposes, if a multimer is requested it's made later
    fasta_creation(monomer_file_name,sequence,1)
    if subunits != 1:
        fasta_creation(multimer_name,sequence,subunits)


def multiple_sequence_alignment(list_of_fastas, fasta_for_alignment, new_alignment_file, reference_protein_fasta='/gpfs/gpfs0/scratch/jws6pq/BridCMfiles/SARS2.fasta'):
    """Creates a multiple sequence alignment using muscle and a concatenated fasta file with a reference fasta as the base,
     joined with all fastas specified in list_of_fastas."""
    from os import system
    # Creating a copy of the fasta file of a reference_protein_fasta to be added into
    system(f'cp  {reference_protein_fasta} {fasta_for_alignment}')
    # Concatenating the fasta files found in list_of_fasta
    for fasta in list_of_fastas:
        system(f"cat {fasta} >> {fasta_for_alignment}")
    # Using muscle to perform the alignment
    system(f'module load gcc/9.2.0 && module load muscle/3.8.31 && muscle -in {fasta_for_alignment} -out {new_alignment_file}')


def alignment_finder(alignment_file, sequence_of_interest,comparison_protein,reference_protein='6vsb_B'):
    """Takes a fasta style alignment and a sequence_of_interest from a reference_protein and returns the sequence of the
    comparison_protein that is outlined in the boundaries of the sequence_of_interest, as well as the python index boundaries
     for the found_alignment of the comparison_protein. reference_protein and comparison_protein must the names following
      '>' in the alignment file"""
    #Must be in fasta format
    alignment = open(alignment_file, "r").read().split('>')
    # Splitting up the sequences names and sequences into a dictionary
    sequence_dictionary = {sequence.split('\n')[0]: ''.join(sequence.split('\n')[1:]) for sequence in alignment if
                           len(sequence) != 0}
    # Recording the specified sequences as variables
    reference_sequence=sequence_dictionary[reference_protein]
    comparison_sequence=sequence_dictionary[comparison_protein]
    # Matching python indexing for the indexing from the alignment with some amount of '-' and indexing in the regular sequence
    reference_sequence_indexing=[ind for ind, x in enumerate(reference_sequence) if x != '-']
    # Creating a regular sequence without '-'
    no_gap_reference_sequence=''.join([x for ind, x in enumerate(reference_sequence) if x != '-'])
    # Boundaries are given in python index
    # Additionally the reference_start is the first residue that is spliced,
    # and reference_end is the first residue that's not spliced
    reference_start=reference_sequence_indexing[no_gap_reference_sequence.find(sequence_of_interest)]
    reference_end=reference_sequence_indexing[no_gap_reference_sequence.find(sequence_of_interest) + len(sequence_of_interest)]
    # Pulling the section of the comparison_sequence that overlaps with the sequence_of_interest
    found_alignment = comparison_sequence[reference_start:reference_end].replace('-','')
    # Recording the indexes of the found_alignment
    splice_start=comparison_sequence.replace('-','').find(found_alignment)
    splice_end=splice_start+len(found_alignment)
    return found_alignment,splice_start,splice_end


def run_emboss_needle(new_emboss_file, sequence_one, sequence_two,needle_directory='/scratch/jws6pq/EMBOSS-6.6.0/emboss/needle'):
    from os import system as syst
    syst(f'{needle_directory} -sprotein -gapopen 10 -gapextend 0.5 -outfile {new_emboss_file} '
         f'-asequence asis:{sequence_one} -bsequence asis:{sequence_two}')
