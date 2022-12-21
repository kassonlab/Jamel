def alignment_finder(alignment_file, sequence_of_interest):
    #Must be in CLW format
    alignment=open(alignment_file, 'r').readlines()
    i=3
    first_sequence=''
    second_sequence=''
    while i<(len(alignment)-1):
        x=alignment[i].split()
        first_sequence+=x[1].rstrip()
        i+=1
        x = alignment[i].split()
        second_sequence+=x[1].rstrip()
        i+=3
    first_seq_indexing=[ind for ind, x in enumerate(first_sequence) if x != '-']
    nogap_first_sequence=''.join([x for ind, x in enumerate(first_sequence) if x != '-'])
    #Splice numbers are given in python index so add 1 for actual residue position
    #Additionally the SpliceStart is the first residue that is spliced,
    # and End is the last residue thats spliced
    SpliceStart=first_seq_indexing[nogap_first_sequence.find(sequence_of_interest)]
    SpliceEnd=first_seq_indexing[nogap_first_sequence.find(sequence_of_interest) + len(sequence_of_interest)]
    return SpliceStart, SpliceEnd
#write if statement to do reverse if there is an error
def RunEmboss(AlignmentFile,Protein,Domain,SequenceofInterest):
    import os
    #Must be in CLW format
    Alignment=open(AlignmentFile,'r').readlines()
    i=3
    FirstSequence=''
    SecondSequence=''
    while i<(len(Alignment)-1):
        x=Alignment[i].split()
        FirstSequence+=x[1].rstrip()
        i+=1
        x = Alignment[i].split()
        SecondSequence+=x[1].rstrip()
        i+=3
    FirstSeqIndexing=[ind for ind, x in enumerate(FirstSequence) if x != '-']
    NogapFirstSequence=''.join([x for ind, x in enumerate(FirstSequence) if x != '-'])
    #Splice numbers are given in python index so add 1 for actual residue position
    #Additionally the SpliceStart is the first residue that is spliced, and End is the first residue of the next domain
    SpliceStart=FirstSeqIndexing[NogapFirstSequence.find(SequenceofInterest)]
    SpliceEnd=FirstSeqIndexing[NogapFirstSequence.find(SequenceofInterest)+len(SequenceofInterest)]
    FoundAlignment=SecondSequence[SpliceStart:SpliceEnd].replace('-','')
    os.system('/scratch/jws6pq/EMBOSS-6.6.0/emboss/needle -sprotein -gapopen 10 -gapextend 0.5 -outfile /scratch/jws6pq/Notebook/Emboss/'+Protein+Domain+'.emboss -asequence asis:'+SequenceofInterest+' -bsequence asis:'+FoundAlignment)
    os.system('/scratch/jws6pq/EMBOSS-6.6.0/emboss/needle -sprotein -gapopen 10 -gapextend 0.5 -outfile /scratch/jws6pq/Notebook/Overall/'+Protein+Domain+'.emboss -asequence asis:'+SequenceofInterest+' -bsequence asis:'+FoundAlignment)
    return SpliceStart, SpliceEnd