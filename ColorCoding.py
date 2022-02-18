#def BlosumColorCoding
import Bio.Align.substitution_matrices
import numpy as np
Blosum62Matrix=Bio.Align.substitution_matrices.load('BLOSUM62')
NumberofSequences=4
#Alignment in FASTA format
Sequences=open('SARS2wEverything.aln',"r")
IDK=''.join([x for x in Sequences if  x[0]!='>'  if x!='']).rstrip().strip().replace('\n','').replace(' ','')
DictionaryofSequences={}
SequenceLength=int(len(IDK)/NumberofSequences)
i=1
j=0
while i<=NumberofSequences:
    DictionaryofSequences[i]=IDK[j:SequenceLength+j]
    j+=SequenceLength
    i+=1
i=1
j=0
ConservationScore=np.empty(SequenceLength)
while i<=NumberofSequences:

print(Blosum62Matrix)
