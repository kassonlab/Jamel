def sequence_splice(fasta_file, boundary_tuple, python_index='Yes'):
    """Takes a fasta sequence and returns the section of the sequence between indexes specified by the boundary one and two,
    as well as the sequence with the specified section replaced with a '-'.
    ABCDEFGH, boundary_one=0, boundary_two=3 Returns ABC and -DEFGH
    The residue at boundary_two is the first residue not included in the splice. If you're using alignment_finder, this is already accounted for."""
    if python_index == 'No':
        boundary_two = boundary_tuple[1]-1
        boundary_one = boundary_tuple[0]-1
    fasta = open(fasta_file, "r")
    sequence = ''.join([x for x in fasta if x[0] != '>' if x != '']).strip().replace('\n', '')
    # The spliced region between the 2 specified boundaries is the first sequence
    # in the list followed by the sequence with the spliced region replace by a '-'
    spliced_region=''.join(sequence[boundary_tuple[0]:boundary_tuple[1]])
    non_spliced_region=sequence.replace(spliced_region,'-')
    return spliced_region,non_spliced_region


def chimera_sequence_creation(section_being_spliced_in, marked_sequence, mark_in_the_sequence='-'):
    """Returns a completed sequence after replacing the '-' in a marked sequence with another sequence fragment"""
    chimera_sequence=marked_sequence.replace(mark_in_the_sequence, section_being_spliced_in)
    return chimera_sequence


def fasta_creation(file_name, sequence, subunits):
    """Creates a fasta file with the given file_name, and replicates the sequence within it the specified number of times
    to create a homo multimer if subunits is greater than 1."""
    from pathlib import Path
    file = open(file_name, 'w')
    for x in range(subunits):
        file.write('>{0}\n{1}\n'.format(Path(file_name).stem,sequence))
    file.close()
