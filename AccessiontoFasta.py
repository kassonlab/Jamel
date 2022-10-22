def AccessionNumbertoFasta(Protein,Accession,Subunits=3):
    from Bio import Entrez
    Entrez.email = 'jws6pq@virginia.edu'
    handle = Entrez.efetch(db='protein', id=Accession,retmode='text',rettype='fasta').readlines()
    Sequence=''.join([x for x in handle if  x[0]!='>'  if x!='']).rstrip().strip().replace('\n','').replace(' ','')
    Monomerfile=open('/gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/'+Protein+'.fasta', 'w')
    Monomerfile.write('>'+Protein+'\n'+Sequence+'\n')
    file = open('/gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/'+str(Subunits)+'mer'+Protein+'.fasta', 'w')
    file.write('>'+Protein+'\n'+Sequence+'\n')
    file.close()
    for x in range(Subunits-1):
        file = open('/gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/'+str(Subunits)+'mer'+Protein+'.fasta', 'a')
        file.write('>'+Protein+'\n'+Sequence+'\n')
        file.close()
def FastatoAlignmentFinder(Protein):
    import os
    os.system('cp  /scratch/jws6pq/CMfiles/SARS2.fasta  /gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/'+Protein+"onSARS2.fasta")
    os.system("cat /gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/"+Protein+".fasta >> /gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/"+Protein+"onSARS2.fasta")
    os.system('module load gcc/9.2.0 && module load muscle/3.8.31 && muscle -in ' + Protein+'onSARS2.fasta -clw -out /scratch/jws6pq/Notebook/Alignment/'+ Protein+'onSARS2.aln')
def MultipleSequenceAlignment(ProteinList):
    import os
    List=[line.split()[-1] for line in open(ProteinList,'r').readlines()]
    ListLength=len(List)
    i=0
    os.system('cp  /gpfs/gpfs0/scratch/jws6pq/BridCMfiles/SARS2.fasta  /gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/SARS2wEverything.fasta')
    while i<ListLength:
        os.system("cat /gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/" + List[i] + ".fasta >> /gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/SARS2wEverything.fasta")
        i+=1
    os.system('module load gcc/9.2.0 && module load muscle/3.8.31 && muscle -in /gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/SARS2wEverything.fasta -out /scratch/jws6pq/Notebook/Alignment/SARS2wEverything.aln')
def pdb2fasta(Protein):
    import os
    from Bio import SeqIO
    os.system('module load gcc/9.2.0 && module load muscle/3.8.31 && muscle -in /gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/SARS2wS1.fasta -out /scratch/jws6pq/Notebook/Alignment/SARS2wS1.fasta.aln')


    PDBFile = 'SARS2w'+Protein+'RBD.pdb'
    with open(PDBFile, 'r') as pdb_file:
        for record in SeqIO.parse(pdb_file, 'pdb-atom'):
            print('>' + record.id)
            print(record.seq)

ProteinList=open('List','r').readlines()
for line in ProteinList:
    ProteinInfo=line.split()
    AccessionNumbertoFasta(ProteinInfo[-1],ProteinInfo[0])
MultipleSequenceAlignment('List')