def AlignmentFinder(AlignmentFile):
    import numpy as np
    import Bio
    Alignment=open(AlignmentFile,'r').readlines()
    i=3
    FirstSequence=''
    SecondSequence=''
    SequenceofInterest='ESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSYNYLYRNLKPFERDISTEIYNCYFPLQSYGFQPTVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLT'
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
    from Bio.Emboss.Applications import WaterCommandline
    SequenceSimilarity=Bio.Emboss.Applications.WaterCommandline(cmd='water',stdout='true',asequence=SequenceofInterest,bsequence=SequenceofInterest,similarity='true',sprotein='true',outfile='test.txt',gapopen=10,gapextend=0.5)
    print(SequenceSimilarity)
    return SpliceStart, SpliceEnd
#write if statement to do reverse if there is an error
AlignmentFinder('QuailHKU30onSARS2.aln')