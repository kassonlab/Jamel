

def DomainExchange(Fastafile1,Fastafile2,Boundary3,Boundary4,Boundary1=224,Boundary2=425,Domain='RBD',Subunits=3):
    import numpy as np
    import os
    Fasta1=open(Fastafile1,"r")
    Fasta2=open(Fastafile2,"r")
    Protein1=os.path.basename(Fasta1).replace('.fasta','')
    Protein2=os.path.basename(Fasta2).replace('.fasta','')

    Title=Protein1+'w'+Protein2+Domain+'.fasta'
    Title2=Protein2+'w'+Protein1+Domain+'.fasta'

    Sequence1=''.join([x for x in Fasta1 if  x[0]!='>'  if x!='']).rstrip().strip().replace('\n','').replace(' ','')
    Sequence2=''.join([x for x in Fasta2 if  x[0]!='>'  if x!='']).rstrip().strip().replace('\n','').replace(' ','')

    Sections=np.empty((3,2), dtype=object)

    Sections[0,0]=''.join(Sequence1[0:Boundary1])
    Sections[1,0]=''.join(Sequence1[Boundary1:Boundary2])
    Sections[2,0]=''.join(Sequence1[Boundary2:])

    Sections[0, 1] = ''.join(Sequence2[0:Boundary3])
    Sections[1, 1] = ''.join(Sequence2[Boundary3:Boundary4])
    Sections[2, 1] = ''.join(Sequence2[Boundary4:])
    
    file = open("/scratch/jws6pq/Notebook/Finished/Fastas/"+str(Subunits)+'mer' +Title, 'w')
    file.write('>'+Protein1 + 'w' + Protein2 + Domain + '\n' + Sections[0, 0] + Sections[1, 1] + Sections[2, 0]+'\n')
    #file2 = open("/scratch/jws6pq/Notebook/Finished/Fastas/" +Title2, 'w')
    #file2.write('>' +Protein2 + 'w' + Protein1 + Domain + '\n' + Sections[0, 1] + Sections[1, 0] + Sections[2, 1]+'\n')
    for x in range(Subunits):
        #file2 = open(str(Subunits)+'mer' +Title, 'a')
        #file2.write('>'+Protein2+'w'+Protein1+Domain+ '\n' + Sections[0,1]+Sections[1,0]+Sections[2,1]+'\n')
        file.close()
        #file2.close()
    AlphaFoldEntry= "/scratch/jws6pq/Notebook/Finished/Fastas/" +str(Subunits)+'mer'+ Fasta2.name + ','
    AlphaFoldEntry2="/scratch/jws6pq/Notebook/Finished/Fastas/"+str(Subunits)+'mer'+Title+','
    return AlphaFoldEntry,AlphaFoldEntry2
