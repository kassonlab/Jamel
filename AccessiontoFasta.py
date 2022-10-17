def AccessionNumbertoFasta(Protein,Accession,Subunits=3):
    from Bio import Entrez
    Entrez.email = 'jws6pq@virginia.edu'
    handle = Entrez.efetch(db='protein', id=Accession,retmode='text',rettype='fasta').readlines()
    Sequence=''.join([x for x in handle if  x[0]!='>'  if x!='']).rstrip().strip().replace('\n','').replace(' ','')
    file = open('/scratch/jws6pq/Notebook/Finished/Fastas/'+str(Subunits)+'mer'+Protein+'.fasta', 'w')
    file.write('>'+Protein+'\n'+Sequence+'\n')
    file.close()
    for x in range(Subunits-1):
        file = open('/scratch/jws6pq/Notebook/Finished/Fastas/'+str(Subunits)+'mer'+Protein+'.fasta', 'a')
        file.write('>'+Protein+'\n'+Sequence+'\n')
        file.close()
def FastatoAlignmentFinder(Protein):
    import os
    os.system('cp  /scratch/jws6pq/CMfiles/SARS2.fasta  /scratch/jws6pq/Notebook/Finished/Fastas/'+Protein+"onSARS2.fasta")
    os.system("cat /scratch/jws6pq/Notebook/Finished/Fastas/"+Protein+".fasta >> /scratch/jws6pq/Notebook/Finished/Fastas/"+Protein+"onSARS2.fasta")
    os.system('module load gcc/9.2.0 && module load muscle/3.8.31 && muscle -in ' + Protein+'onSARS2.fasta -clw -out /scratch/jws6pq/Notebook/Alignment/'+ Protein+'onSARS2.aln')
def MultipleSequenceAlignment(ProteinList):
    import os
    List=[line.split()[-1] for line in open(ProteinList,'r').readlines()]
    ListLength=len(List)
    i=0
    os.system('cp  /scratch/jws6pq/CMfiles/SARS2.fasta  /scratch/jws6pq/Notebook/Finished/Fastas/SARS2wS1.fasta')
    while i<ListLength:
        os.system("cat /scratch/jws6pq/Notebook/Finished/Fastas/SARS2w" + List[i] + "RBD.fasta >> /scratch/jws6pq/Notebook/Finished/Fastas/SARS2wRBD.fasta")
        i+=1
    os.system('module load gcc/9.2.0 && module load muscle/3.8.31 && muscle -in /scratch/jws6pq/Notebook/Finished/Fastas/SARS2wRBD.fasta -out /scratch/jws6pq/Notebook/Alignment/SARS2wRBD.fasta.aln')
def pdb2fasta(Protein):
    import sys
    from Bio import SeqIO
    os.system('module load gcc/9.2.0 && module load muscle/3.8.31 && muscle -in /scratch/jws6pq/Notebook/Finished/Fastas/SARS2wS1.fasta -out /scratch/jws6pq/Notebook/Alignment/SARS2wS1.fasta.aln')


    PDBFile = 'SARS2w'+Protein+'RBD.pdb'
    with open(PDBFile, 'r') as pdb_file:
        for record in SeqIO.parse(pdb_file, 'pdb-atom'):
            print('>' + record.id)
            print(record.seq)
# MultipleSequenceAlignment('List')
pdb2fasta('Avian')