def AlignmentFinder(AlignmentFile,SequenceofInterest):
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
    #Additionally the SpliceStart is the first residue that is spliced, and End is the last residue thats spliced
    SpliceStart=FirstSeqIndexing[NogapFirstSequence.find(SequenceofInterest)]
    SpliceEnd=FirstSeqIndexing[NogapFirstSequence.find(SequenceofInterest)+len(SequenceofInterest)]
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
    #Additionally the SpliceStart is the first residue that is spliced, and End is the last residue thats spliced
    SpliceStart=FirstSeqIndexing[NogapFirstSequence.find(SequenceofInterest)]
    SpliceEnd=FirstSeqIndexing[NogapFirstSequence.find(SequenceofInterest)+len(SequenceofInterest)]
    FoundAlignment=SecondSequence[SpliceStart:SpliceEnd].replace('-','')
    os.system('/sfs/lustre/bahamut/scratch/jws6pq/EMBOSS-6.6.0/emboss/needle -sprotein -gapopen 10 -gapextend 0.5 -outfile /sfs/lustre/bahamut/scratch/jws6pq/Notebook/Emboss/'+Protein+Domain+'.emboss -asequence asis:'+SequenceofInterest+' -bsequence asis:'+FoundAlignment)
    os.system('/sfs/lustre/bahamut/scratch/jws6pq/EMBOSS-6.6.0/emboss/needle -sprotein -gapopen 10 -gapextend 0.5 -outfile /sfs/lustre/bahamut/scratch/jws6pq/Notebook/Overall/'+Protein+Domain+'.emboss -asequence asis:'+SequenceofInterest+' -bsequence asis:'+FoundAlignment)
    return SpliceStart, SpliceEnd