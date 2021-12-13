def AlignmentFinder(AlignmentFile,Protein):
    import numpy as np
    import Bio
    import os
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
    os.system('/sfs/lustre/bahamut/scratch/jws6pq/EMBOSS-6.6.0/emboss/needle -sprotein -gapopen 10 -gapextend 0.5 -outfile /sfs/lustre/bahamut/scratch/jws6pq/Notebook/Emboss/'+Protein+'.emboss -asequence asis:'+SequenceofInterest+' -bsequence asis:'+FoundAlignment)
    return SpliceStart, SpliceEnd
#write if statement to do reverse if there is an error
