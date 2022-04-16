# def BlosumColorCoding():
import numpy as np
from pymol import cmd
Blosum62Matrix = np.array([['','A','R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W','Y', 'V', 'B', 'Z', 'X', '-'],
                       ['A', 4.0, -1.0, -2.0, -2.0, 0.0, -1.0, -1.0, 0.0, -2.0, -1.0, -1.0, -1.0, -1.0, -2.0, -1.0, 1.0,0.0, -3.0, -2.0, 0.0, -2.0, -1.0, 0.0, -4.0],
                       ['R', -1.0, 5.0, 0.0, -2.0, -3.0, 1.0, 0.0, -2.0, 0.0, -3.0, -2.0, 2.0, -1.0, -3.0, -2.0, -1.0,-1.0, -3.0, -2.0, -3.0, -1.0, 0.0, -1.0, -4.0],
                       ['N', -2.0, 0.0, 6.0, 1.0, -3.0, 0.0, 0.0, 0.0, 1.0, -3.0, -3.0, 0.0, -2.0, -3.0, -2.0, 1.0, 0.0,-4.0, -2.0, -3.0, 3.0, 0.0, -1.0, -4.0],
                       ['D', -2.0, -2.0, 1.0, 6.0, -3.0, 0.0, 2.0, -1.0, -1.0, -3.0, -4.0, -1.0, -3.0, -3.0, -1.0, 0.0,-1.0, -4.0, -3.0, -3.0, 4.0, 1.0, -1.0, -4.0],
                       ['C', 0.0, -3.0, -3.0, -3.0, 9.0, -3.0, -4.0, -3.0, -3.0, -1.0, -1.0, -3.0, -1.0, -2.0, -3.0,-1.0, -1.0, -2.0, -2.0, -1.0, -3.0, -3.0, -2.0, -4.0],
                       ['Q', -1.0, 1.0, 0.0, 0.0, -3.0, 5.0, 2.0, -2.0, 0.0, -3.0, -2.0, 1.0, 0.0, -3.0, -1.0, 0.0,-1.0, -2.0, -1.0, -2.0, 0.0, 3.0, -1.0, -4.0],
                       ['E', -1.0, 0.0, 0.0, 2.0, -4.0, 2.0, 5.0, -2.0, 0.0, -3.0, -3.0, 1.0, -2.0, -3.0, -1.0, 0.0,-1.0, -3.0, -2.0, -2.0, 1.0, 4.0, -1.0, -4.0],
                       ['G', 0.0, -2.0, 0.0, -1.0, -3.0, -2.0, -2.0, 6.0, -2.0, -4.0, -4.0, -2.0, -3.0, -3.0, -2.0, 0.0,-2.0, -2.0, -3.0, -3.0, -1.0, -2.0, -1.0, -4.0],
                       ['H', -2.0, 0.0, 1.0, -1.0, -3.0, 0.0, 0.0, -2.0, 8.0, -3.0, -3.0, -1.0, -2.0, -1.0, -2.0, -1.0,-2.0, -2.0, 2.0, -3.0, 0.0, 0.0, -1.0, -4.0],
                       ['I', -1.0, -3.0, -3.0, -3.0, -1.0, -3.0, -3.0, -4.0, -3.0, 4.0, 2.0, -3.0, 1.0, 0.0, -3.0, -2.0,-1.0, -3.0, -1.0, 3.0, -3.0, -3.0, -1.0, -4.0],
                       ['L', -1.0, -2.0, -3.0, -4.0, -1.0, -2.0, -3.0, -4.0, -3.0, 2.0, 4.0, -2.0, 2.0, 0.0, -3.0, -2.0,-1.0, -2.0, -1.0, 1.0, -4.0, -3.0, -1.0, -4.0],
                       ['K', -1.0, 2.0, 0.0, -1.0, -3.0, 1.0, 1.0, -2.0, -1.0, -3.0, -2.0, 5.0, -1.0, -3.0, -1.0, 0.0,-1.0, -3.0, -2.0, -2.0, 0.0, 1.0, -1.0, -4.0],
                       ['M', -1.0, -1.0, -2.0, -3.0, -1.0, 0.0, -2.0, -3.0, -2.0, 1.0, 2.0, -1.0, 5.0, 0.0, -2.0, -1.0,-1.0, -1.0, -1.0, 1.0, -3.0, -1.0, -1.0, -4.0],
                       ['F', -2.0, -3.0, -3.0, -3.0, -2.0, -3.0, -3.0, -3.0, -1.0, 0.0, 0.0, -3.0, 0.0, 6.0, -4.0, -2.0,-2.0, 1.0, 3.0, -1.0, -3.0, -3.0, -1.0, -4.0],
                       ['P', -1.0, -2.0, -2.0, -1.0, -3.0, -1.0, -1.0, -2.0, -2.0, -3.0, -3.0, -1.0, -2.0, -4.0, 7.0,-1.0, -1.0, -4.0, -3.0, -2.0, -2.0, -1.0, -2.0, -4.0],
                       ['S', 1.0, -1.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0, -2.0, -2.0, 0.0, -1.0, -2.0, -1.0, 4.0,1.0, -3.0, -2.0, -2.0, 0.0, 0.0, 0.0, -4.0],
                       ['T', 0.0, -1.0, 0.0, -1.0, -1.0, -1.0, -1.0, -2.0, -2.0, -1.0, -1.0, -1.0, -1.0, -2.0, -1.0,1.0, 5.0, -2.0, -2.0, 0.0, -1.0, -1.0, 0.0, -4.0],
                       ['W', -3.0, -3.0, -4.0, -4.0, -2.0, -2.0, -3.0, -2.0, -2.0, -3.0, -2.0, -3.0, -1.0, 1.0, -4.0,-3.0, -2.0, 11.0, 2.0, -3.0, -4.0, -3.0, -2.0, -4.0],
                       ['Y', -2.0, -2.0, -2.0, -3.0, -2.0, -1.0, -2.0, -3.0, 2.0, -1.0, -1.0, -2.0, -1.0, 3.0, -3.0,-2.0, -2.0, 2.0, 7.0, -1.0, -3.0, -2.0, -1.0, -4.0],
                       ['V', 0.0, -3.0, -3.0, -3.0, -1.0, -2.0, -2.0, -3.0, -3.0, 3.0, 1.0, -2.0, 1.0, -1.0, -2.0, -2.0,0.0, -3.0, -1.0, 4.0, -3.0, -2.0, -1.0, -4.0],
                       ['B', -2.0, -1.0, 3.0, 4.0, -3.0, 0.0, 1.0, -1.0, 0.0, -3.0, -4.0, 0.0, -3.0, -3.0, -2.0, 0.0,-1.0, -4.0, -3.0, -3.0, 4.0, 1.0, -1.0, -4.0],
                       ['Z', -1.0, 0.0, 0.0, 1.0, -3.0, 3.0, 4.0, -2.0, 0.0, -3.0, -3.0, 1.0, -1.0, -3.0, -1.0, 0.0,-1.0, -3.0, -2.0, -2.0, 1.0, 4.0, -1.0, -4.0],
                       ['X', 0.0, -1.0, -1.0, -1.0, -2.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -2.0,0.0, 0.0, -2.0, -1.0, -1.0, -1.0, -1.0, -1.0, -4.0],
                       ['-', -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0,-4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0, 1.0],
                       ], dtype=object)
# Alignment in FASTA format. Make sure your benchmark sequence is first
Sequences =open('SARS2stable.aln', "r").readlines()
SequenceCount=''.join(Sequences)
NumberofSequences = SequenceCount.count('>')

#This is taking all the sequences of the fasta file and assigning them to a dictionary
EntireFasta= ''.join([x for x in Sequences if x[0] != '>']).rstrip().strip().replace('\n', '').replace(' ', '')
DictionaryofSequences = {}
SequenceLength = int(len(EntireFasta) / NumberofSequences)
i = 1
j = 0
while i <= NumberofSequences:
    DictionaryofSequences[i] = EntireFasta[j:SequenceLength + j]
    j += SequenceLength
    i += 1

#This is comparing the each residue of SARS2 and comparing it to the aligned sequence residues. Using BLOSUM matrix its averaging the substition cost for each residue
Residueindex = 0
Sequenceindex = 2
ConservationScore = np.zeros((3,SequenceLength), dtype=object)

while Residueindex<SequenceLength:
    SARS2Residue=np.where(Blosum62Matrix[:,0]==DictionaryofSequences[1][Residueindex])[0]
    ConservationScore[0,Residueindex]=DictionaryofSequences[1][Residueindex]
    while Sequenceindex<=NumberofSequences:
        ComparisonResidue=np.where(Blosum62Matrix[0,:]==DictionaryofSequences[Sequenceindex][Residueindex])[0]
        ConservationScore[1,Residueindex]+=int(Blosum62Matrix[SARS2Residue,ComparisonResidue])
        Sequenceindex+=1
    ConservationScore[1,Residueindex]=ConservationScore[1,Residueindex]/NumberofSequences
    Residueindex+=1
    Sequenceindex=2

#This is removing all the colums with dashes that were added in the alignment
TruncatedConservationScore=np.delete(ConservationScore,np.where(ConservationScore[0]=='-'),axis=1)
# def ColorCoding():

#This is calculating how well conserved each SARS residue is based on how far it is away from the optimal substitution which is substitution of the same residue. The fraction is assigned to the residue and assigned a color based on the spectrum provided.
Residueindex=0
for x in TruncatedConservationScore[1]:
    MatrixRow=Blosum62Matrix[np.where(Blosum62Matrix[:, 0] == TruncatedConservationScore[0,Residueindex])[0][0], 1:]
    TruncatedConservationScore[2,Residueindex]=1-((int(max(MatrixRow)))-x)/np.ptp(MatrixRow)
    cmd.alter('resi '+str(Residueindex+1),'b='+str(TruncatedConservationScore[2,Residueindex]))
    Residueindex+=1
#Lower values are assigned the first color listed, higher values latter color
cmd.spectrum('b','blue_white_green','6VSB_B')
