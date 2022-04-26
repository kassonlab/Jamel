def BlosumColorCoding(Alignmentfile):
    import numpy as np
    from collections import Counter
    import math
    import sys
    # from pymol import cmd
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
    Sequences =open(Alignmentfile, "r").readlines()
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
    #Storing the top sequence of the alignment as the reference sequence
    SARS2Sequence=DictionaryofSequences[1]
    Residueindex=0
    Sequenceindex = 2
    SubstitutionScoresperResidue = list(np.ones((NumberofSequences-1), dtype=int))
    ShannonsEntropy=np.zeros((2,SequenceLength), dtype=object)
    while Residueindex<SequenceLength:
        if SARS2Sequence[Residueindex]=='-':
            ShannonsEntropy[0, Residueindex] = SARS2Sequence[Residueindex]
            Residueindex += 1
        else:
            SARS2Residueindex=np.where(Blosum62Matrix[:,0]==SARS2Sequence[Residueindex])[0]
            ShannonsEntropy[0,Residueindex]=SARS2Sequence[Residueindex]
            while Sequenceindex<=NumberofSequences:
                ComparisonResidueindex=np.where(Blosum62Matrix[0,:]==DictionaryofSequences[Sequenceindex][Residueindex])[0]
                SubstitutionScoresperResidue[Sequenceindex-2]=int(Blosum62Matrix[SARS2Residueindex,ComparisonResidueindex])
                Sequenceindex+=1
            MatrixValuesHistogram=Counter(SubstitutionScoresperResidue)
            NumberofBins=len(MatrixValuesHistogram)

            for key, value in MatrixValuesHistogram.items():
                if NumberofBins==1:
                    ShannonsEntropy[1, Residueindex]=0
                else:
                    probability=value/(NumberofSequences-1)
                    ShannonsEntropy[1,Residueindex]+=probability*math.log(probability,NumberofBins)
            ShannonsEntropy[1, Residueindex]=(ShannonsEntropy[1,Residueindex])/-1
            Residueindex+=1
            Sequenceindex=2
    np.set_printoptions(threshold=sys.maxsize)
    ShannonsEntropy=np.delete(ShannonsEntropy,np.where(ShannonsEntropy[0]=='-'),axis=1)
    print(ShannonsEntropy)
    # for x in TruncatedConservationScore[1]:
    #     cmd.alter('resi '+str(Residueindex+1),'b='+str(TruncatedConservationScore[2,Residueindex]))
    #     Residueindex+=1
    # #Lower values are assigned the first color listed, higher values latter color
    # cmd.spectrum('b','blue_green','6VSB_B')
BlosumColorCoding('Random Sequence')