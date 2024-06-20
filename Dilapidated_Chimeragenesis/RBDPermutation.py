import numpy as np
import numpy.random.mtrand

file=open('PANG.fasta',"r")
Fasta1=file.readlines()
for line in Fasta1:
    Sequence1=line
Seqlength=len(Sequence1)
Sections=np.empty((3), dtype=object)

Sections[0]=''
Sections[1]=''
Sections[2]=''
i=0
j=0
for letter in Sequence1:
    i+=1
    if i<=Seqlength*(331/1141):
        Sections[0]+=str(letter)
    elif Seqlength*(331/1141)<i<=Seqlength*(524/1141):
        Sections[1] += str(letter)
    elif Seqlength*(524/1141)<i:
        Sections[2] += str(letter)
RBDLength=len(Sections[1])
RBDSequence=np.empty(RBDLength, dtype=object)
i=0
for letter in Sections[1]:
    RBDSequence[i]=letter
    i+=1

Randomized=numpy.random.mtrand.permutation(RBDSequence)

NewSequence=np.empty((2), dtype=object)
NewSequence[0]='>Chimera'
NewSequence[1]=Sections[0]
i=0
for i in range(len(Randomized)):
    NewSequence[1]+=Randomized[i]
    i+=1
NewSequence[1]+=Sections[2]
print(NewSequence)
print(Randomized)

np.savetxt('PANGRBDPermute.fasta',NewSequence,fmt="%s",delimiter=" ")
exit()