def AlignmentFinder(AlignmentFile,Protein):
    import os
    Alignment=open(AlignmentFile,'r').readlines()
    i=3
    FirstSequence=''
    SecondSequence=''
    SequenceofInterest='TSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSYNYLYRNLKPFERDISTEIYNCYFPLQSYGFQPTVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVN'
    while i<(len(Alignment)-1):
        x=Alignment[i].split()
        FirstSequence+=x[1].rstrip()
        i+=1
        x = Alignment[i].split()
        SecondSequence+=x[1].rstrip()
        i+=3
    FirstSeqIndexing=[ind for ind, x in enumerate(FirstSequence) if x != '-']
    NogapFirstSequence=''.join([x for ind, x in enumerate(FirstSequence) if x != '-'])
    SpliceStart=FirstSeqIndexing[NogapFirstSequence.find(SequenceofInterest)]
    SpliceEnd=FirstSeqIndexing[NogapFirstSequence.find(SequenceofInterest)+len(SequenceofInterest)]
    FoundAlignment=SecondSequence[SpliceStart:SpliceEnd].replace('-','')
    SpliceStart=SecondSequence.replace('-','').find(FoundAlignment)
    SpliceEnd=len(FoundAlignment)+SpliceStart
    os.system('/sfs/lustre/bahamut/scratch/jws6pq/EMBOSS-6.6.0/emboss/needle -sprotein -gapopen 10 -gapextend 0.5 -outfile /sfs/lustre/bahamut/scratch/jws6pq/Notebook/Emboss/'+Protein+'S1.emboss -asequence asis:'+SequenceofInterest+' -bsequence asis:'+FoundAlignment)
    os.system('/sfs/lustre/bahamut/scratch/jws6pq/EMBOSS-6.6.0/emboss/needle -sprotein -gapopen 10 -gapextend 0.5 -outfile /sfs/lustre/bahamut/scratch/jws6pq/Notebook/Overall/'+Protein+'S1.emboss -asequence asis:'+SequenceofInterest+' -bsequence asis:'+FoundAlignment)
    return SpliceStart, SpliceEnd
AlignmentFinder('FullSARSonSARS2.aln','SARS')