def AlignmentFinder(AlignmentFile,SequenceofInterest):
    import numpy as np
    Alignment=open('SARS2onFullMERS.aln','r').readlines()
    i=3
    FirstSequence=''
    SecondSequence=''
    SequenceofInterest='ESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSYNYLYRNLKPFERDISTEIYNCYFPLQSYGFQPTVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLT'
    while i<(len(Alignment)-1):
        x=Alignment[i].split()
        print(x)
        FirstSequence+=x[1].rstrip()
        i+=1
        x = Alignment[i].split()
        print(x)
        SecondSequence+=x[1].rstrip()
        i+=3
    print()
    FirstSeqIndexing=[ind for ind, x in enumerate(FirstSequence) if x != '-']
    print(FirstSeqIndexing)
    NogapFirstSequence=''.join([x for ind, x in enumerate(FirstSequence) if x != '-'])
    print(NogapFirstSequence.find(SequenceofInterest))
    SpliceStart=FirstSeqIndexing[NogapFirstSequence.find(SequenceofInterest)]
    print(SpliceStart)
    SpliceEnd=FirstSeqIndexing[NogapFirstSequence.find(SequenceofInterest)+len(SequenceofInterest)]
    print(SpliceEnd)
    FoundAlignment=SecondSequence[SpliceStart:SpliceEnd].replace('-','')
    print(FoundAlignment)
#write if statement to do reverse if there is an error
