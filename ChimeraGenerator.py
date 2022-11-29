def SequenceSplice(Fastafile, Boundary1, Boundary2,Pythonindexforboundary='Yes'):
    import os

    if Pythonindexforboundary=='No':
        Boundary1+=-1
        Boundary2+=-1

    Fasta = open(Fastafile, "r")

    Sequence = ''.join([x for x in Fasta if x[0] != '>' if x != '']).rstrip().strip().replace('\n', '').replace(' ','')
    #The spliced region between the 2 specified boundaries is the first sequence in the list followed by the sequence with the spliced region replace by a '-'
    SplicedRegion=''.join(Sequence[Boundary1:Boundary2])
    Nonspliced=Sequence.replace(SplicedRegion,'-')
    return SplicedRegion,Nonspliced
def ChimeraSequenceCreation(SectiontobeSplicedin,MarkedSequence):
    ChimeraSequence=MarkedSequence.replace('-',SectiontobeSplicedin)

    ChimeraSplicetuple=(ChimeraSequence.find('SectiontobeSplicedin'),ChimeraSequence.find('SectiontobeSplicedin')+len(SectiontobeSplicedin))
    Nonsplicetuple=[(ChimeraSequence.find(x),ChimeraSequence.find(x)+len(x)) for x in MarkedSequence.split('-')]
    return ChimeraSequence,ChimeraSplicetuple,Nonsplicetuple


def FastaCreation(Filename,Sequence,Subunits=1):
    file = open(Filename, 'w')
    if '/' in Filename:
        Filename=Filename.split('/')[-1]
    file.write('>' +Filename.replace('.fasta','') + '\n' + Sequence + '\n')
    for x in range(Subunits-1):
        file.write('>' +Filename.replace('.fasta','') + '\n' + Sequence + '\n')
    file.close()
