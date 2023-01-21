def alignment_finder(alignment_file, sequence_of_interest,comparison_protein,reference_protein='6vsb_B'):
    #Must be in fasta format
    alignment = open(alignment_file, "r").read().split('>')
    sequence_dictionary = {sequence.split('\n')[0]: ''.join(sequence.split('\n')[1:]) for sequence in alignment if
                           len(sequence) != 0}
    reference_sequence=sequence_dictionary[reference_protein]
    comparison_sequence=sequence_dictionary[comparison_protein]
    reference_sequence_indexing=[ind for ind, x in enumerate(reference_sequence) if x != '-']
    no_gap_reference_sequence=''.join([x for ind, x in enumerate(reference_sequence) if x != '-'])
    # Boundaries are given in python index
    # Additionally the reference_start is the first residue that is spliced,
    # and reference_end is the first residue that's not spliced
    reference_start=reference_sequence_indexing[no_gap_reference_sequence.find(sequence_of_interest)]
    reference_end=reference_sequence_indexing[no_gap_reference_sequence.find(sequence_of_interest) + len(sequence_of_interest)]
    found_alignment = comparison_sequence[reference_start:reference_end].replace('-','')
    splice_start=comparison_sequence.replace('-','').find(found_alignment)
    splice_end=splice_start+len(found_alignment)
    return found_alignment,splice_start,splice_end


def run_emboss(new_emboss_file, sequence_one, sequence_two):
    from os import system as syst
    syst(f'/scratch/jws6pq/EMBOSS-6.6.0/emboss/needle -sprotein -gapopen 10 -gapextend 0.5 -outfile {new_emboss_file} -asequence asis:{sequence_one} -bsequence asis:{sequence_two}')