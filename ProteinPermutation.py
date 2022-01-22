def Permutation(Protein):
    import numpy as np
    import numpy.random.mtrand

    Fasta=open(Protein+'.fasta',"r").readlines()
    Sequence=list(''.join([x for x in Fasta if  x[0]!='>'  if x!='']).rstrip().strip().replace('\n','').replace(' ',''))
    RandomizedSequence=''.join(numpy.random.mtrand.permutation(Sequence))
    NewSequence=np.empty((2), dtype=object)
    NewSequence[0]='>Chimera'
    NewSequence[1]=RandomizedSequence
    print(NewSequence[1])
    np.savetxt('Permuted'+Protein+'.fasta',NewSequence,fmt="%s",delimiter=" ")
Permutation('SARS2')