def sequence_splice(fasta_file, boundary_one, boundary_two, python_index='Yes'):
    if python_index == 'Yes':
        boundary_two += 1
    else:
        boundary_one += -1

    fasta = open(fasta_file, "r")
    sequence = ''.join([x for x in fasta if x[0] != '>' if x != '']).rstrip().strip().replace('\n', '').replace(' ','')
    # The spliced region between the 2 specified boundaries is the first sequence
    # in the list followed by the sequence with the spliced region replace by a '-'
    spliced_region=''.join(sequence[boundary_one:boundary_two])
    nonspliced=sequence.replace(spliced_region,'-')
    return spliced_region,nonspliced


def chimera_sequence_creation(section_being_spliced_in, marked_sequence, mark_in_the_sequence='-'):
    chimera_sequence=marked_sequence.replace(mark_in_the_sequence, section_being_spliced_in)
    return chimera_sequence


def fasta_creation(file_name, sequence, subunits=3):
    file = open(file_name, 'w')
    if '/' in file_name:
        file_name=file_name.split('/')[-1]
    file.write('>' + file_name.replace('.fasta', '') + '\n' + sequence + '\n')
    for x in range(subunits):
        file.write('>' + file_name.replace('.fasta', '') + '\n' + sequence + '\n')
    file.close()
