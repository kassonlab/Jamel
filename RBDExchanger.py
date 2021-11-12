import numpy as np
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
    NogapSecondSequence = ''.join([x for ind, x in enumerate(SecondSequence) if x != '-'])
    print(NogapFirstSequence.find(SequenceofInterest))
    SpliceStart=FirstSeqIndexing[NogapFirstSequence.find(SequenceofInterest)]
    print(SpliceStart)
    SpliceEnd=FirstSeqIndexing[NogapFirstSequence.find(SequenceofInterest)+len(SequenceofInterest)]
    print(SpliceEnd)
    FoundAlignment=SecondSequence[SpliceStart:SpliceEnd].replace('-','')
    print(FoundAlignment)
def RBDExchange(Fastafile1,Fastafile2,Boundary1,Boundary2):
    file=open(Fastafile1,"r")
    file2=open(Fastafile2,"r")
    Fasta1=file.readlines()
    Fasta2=file2.readlines()
    Protein1=file.name.replace('.fasta','')
    Protein2=file2.name.replace('.fasta','')
    Title=Protein1+'w'+Protein2+'RBD.fasta'
    Path='E:\Fastas\\'  +Title
    Title2=Protein2+'w'+Protein1+'RBD.fasta'
    Path2='E:\Fastas\\'  +Title2
    #Default is 319
    FirstSplice=Boundary1
    #Sefault is 541
    SecondSplice=Boundary2
    Sequence1=''
    Sequence2=''
    for i in range(len(Fasta1)-1):
        Sequence1+=Fasta1[i+1]
    for i in range(len(Fasta2)-1):
        Sequence2+=Fasta2[i+1]
    Seqlength=len(Sequence1)
    Seqlength2=len(Sequence2)
    Sections=np.empty((3,2), dtype=object)
    ProteinList=np.empty((4,3), dtype=object)
    Sections[0,0]=''
    Sections[1,0]=''
    Sections[2,0]=''
    i=0
    j=0
    for letter in Sequence1:
        i+=1
        if i<Seqlength*(FirstSplice/1273):
            Sections[0,0]+=str(letter)
        elif Seqlength*(FirstSplice/1273)<i<=Seqlength*(SecondSplice/1273):
            Sections[1,0] += str(letter)
        elif Seqlength*(SecondSplice/1273)<i:
            Sections[2,0] += str(letter)
    i=0
    for i in range(3):
        Sections[i,1]=''
    i=0
    j=1
    for letter in Sequence2:
        i+=1
        if i<Seqlength2*(FirstSplice/1273):
            Sections[0,1]+=str(letter)
        elif Seqlength2*(FirstSplice/1273)<i<=Seqlength2*(SecondSplice/1273):
            Sections[1,1] += str(letter)
        elif Seqlength2*(SecondSplice/1273)<i:
            Sections[2,1] += str(letter)
    NewSequence=np.empty((2,1), dtype=object)
    NewSequence[0,0]='>Chimera'
    NewSequence[1,0]=Sections[0,0]+Sections[1,1]+Sections[2,0]
    m=0

    np.savetxt(Path,NewSequence,fmt="%s",delimiter=" ")
    NewSequence2=np.empty((2,1), dtype=object)
    NewSequence2[0,0]='>Chimera'
    NewSequence2[1,0]=Sections[0,1]+Sections[1,0]+Sections[2,1]
    np.savetxt(Path2,NewSequence2,fmt="%s",delimiter=" ")
    ProteinList[m,0]=Protein1 + '.plddt'
    ProteinList[m+1,0] =Protein2 + '.plddt'
    ProteinList[m+2,0] = 'E:\Plddt\\' +Protein1 +'w'+Protein2+ '_'+'TruncS1.plddt'
    ProteinList[m+3,0] = 'E:\Plddt\\' +Protein2+'w'+Protein1+ '_'+'TruncS1.plddt'
    ProteinList[m, 1]  = ProteinList[m+2,1]=len(Sections[0,0])
    ProteinList[m+1, 1] = ProteinList[m + 3, 1] =len(Sections[0, 1])
    ProteinList[m, 2] =ProteinList[m+3,2]=len(Sections[1, 0])
    ProteinList[m + 2, 2] =ProteinList[m+1,2]=len(Sections[1, 1])
    np.savetxt('EverythingProteinList.tsv', ProteinList, fmt="%s", delimiter=" ")
    AlphaFoldEntry=''
    AlphaFoldEntry+= "/scratch/jws6pq/Notebook/Finished/" + Title +',' + "/scratch/jws6pq/Notebook/Finished/" + Fastafile2
    return AlphaFoldEntry