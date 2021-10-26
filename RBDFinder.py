#def AlignmentFinder(AlignmentFile):
import numpy as np
Alignment=np.loadtxt('SARS2onFullMERS.aln',dtype=str,skiprows=3,usecols=(1,))
i=0
FirstSequence=''
SecondSequence=''
SequenceofInterest='ESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLT'
while i<(len(Alignment)-1):
    FirstSequence+=Alignment[i].rstrip()
    i+=1
    SecondSequence+=Alignment[i].rstrip()
    i+=2
print(enumerate(SequenceofInterest))
SpliceStart=FirstSequence.index(SequenceofInterest[0:6])
SpliceEnd=FirstSequence.index((SequenceofInterest[-7:-1]))
SecondSequence[SpliceStart:SpliceEnd]
