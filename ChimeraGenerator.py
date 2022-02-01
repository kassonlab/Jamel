import numpy as np

def DomainExchange(Fastafile1,Fastafile2,Boundary3,Boundary4,Boundary1=224,Boundary2=424,Domain='RBD'):
    Fasta1=open(Fastafile1,"r")
    Fasta2=open(Fastafile2,"r")
    Protein1=Fasta1.name.replace('.fasta','')
    Protein2=Fasta2.name.replace('.fasta','')

    Title=Protein1+'w'+Protein2+Domain+'.fasta'
    Title2=Protein2+'w'+Protein1+Domain+'.fasta'

    Sequence1=''.join([x for x in Fasta1 if  x[0]!='>'  if x!='']).rstrip().strip().replace('\n','').replace(' ','')
    Sequence2=''.join([x for x in Fasta2 if  x[0]!='>'  if x!='']).rstrip().strip().replace('\n','').replace(' ','')

    Sections=np.empty((3,2), dtype=object)

    Sections[0,0]=''.join(Sequence1[0:Boundary1-1])
    Sections[1,0]=''.join(Sequence1[Boundary1-1:Boundary2-1])
    Sections[2,0]=''.join(Sequence1[Boundary2-1:])

    Sections[0, 1] = ''.join(Sequence2[0:Boundary3 - 1])
    Sections[1, 1] = ''.join(Sequence2[Boundary3-1:Boundary4 - 1])
    Sections[2, 1] = ''.join(Sequence2[Boundary4-1:])

    NewSequence=np.empty((2,1), dtype=object)
    NewSequence[0,0]='>Chimera'
    NewSequence[1,0]=Sections[0,0]+Sections[1,1]+Sections[2,0]
    np.savetxt(Title,NewSequence,fmt="%s",delimiter=" ")

    NewSequence2=np.empty((2,1), dtype=object)
    NewSequence2[0,0]='>Chimera'
    NewSequence2[1,0]=Sections[0,1]+Sections[1,0]+Sections[2,1]
    #np.savetxt(Title2,NewSequence2,fmt="%s",delimiter=" ")
    AlphaFoldEntry= "/scratch/jws6pq/Notebook/Finished/Fastas/" + Fasta2.name + ','
    AlphaFoldEntry2="/scratch/jws6pq/Notebook/Finished/Fastas/"+Title+','
    return AlphaFoldEntry,AlphaFoldEntry2
