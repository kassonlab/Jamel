def SequenceSplice(Fastafile1, Boundary1, Boundary2,Pythonindexforboundary='Yes'):
    import numpy as np
    import os

    if Pythonindexforboundary=='No':
        Boundary1+=-1
        Boundary2+=-1

    Fasta1 = open(Fastafile1, "r")
    Protein1 = os.path.basename(Fastafile1).replace('.fasta', '')

    Sequence1 = ''.join([x for x in Fasta1 if x[0] != '>' if x != '']).rstrip().strip().replace('\n', '').replace(' ','')

    Sections = np.empty((4, 1), dtype=object)

    Sections[0, 0] = Protein1
    #The spliced region between the 2 specified boundaries is the first sequence in the list followed by what was before and then after the splice
    Sections[1, 0] = ''.join(Sequence1[Boundary1:Boundary2])
    Sections[2, 0] = ''.join(Sequence1[0:Boundary1])
    Sections[3, 0] = ''.join(Sequence1[Boundary2:])

    return [x for x in Sections[:,0] if len(x)>0]

def FastaCreation(Filename,Sequence,Subunits=1):
    file = open(Filename, 'w')
    if '/' in Filename:
        Filename=Filename.split('/')[-1]
    file.write('>' +Filename.replace('.fasta','') + '\n' + Sequence + '\n')
    for x in range(Subunits-1):
        file.write('>' +Filename.replace('.fasta','') + '\n' + Sequence + '\n')
    file.close()
