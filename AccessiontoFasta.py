def AccessionNumbertoFasta(Protein,Accession):
    import numpy as np
    from Bio import Entrez
    Entrez.email = 'jws6pq@virginia.edu'
    handle = Entrez.efetch(db='protein', id=Accession,retmode='text',rettype='fasta').readlines()
    handle= [x for x in handle if x!='']
    np.savetxt('/sfs/lustre/bahamut/scratch/jws6pq/Notebook/Finished/Fastas/'+Protein+'.fasta', handle , fmt="%s", delimiter="")
def FastatoAlignmentFinder(Protein):
    import os
    os.system('cp  /sfs/lustre/bahamut/scratch/jws6pq/CMfiles/SARS2.fasta  /sfs/lustre/bahamut/scratch/jws6pq/Notebook/Finished/Fastas/'+Protein+"onSARS2.fasta")
    os.system("cat /sfs/lustre/bahamut/scratch/jws6pq/Notebook/Finished/Fastas/"+Protein+".fasta >> /sfs/lustre/bahamut/scratch/jws6pq/Notebook/Finished/Fastas/"+Protein+"onSARS2.fasta")
    os.system('module load gcc/9.2.0 && module load muscle/3.8.31 && muscle -in ' + Protein+'onSARS2.fasta -clw -out /sfs/lustre/bahamut/scratch/jws6pq/Notebook/Alignment/'+ Protein+'onSARS2.aln')
def MultipleSequenceAlignment(ProteinList):
    import os
    List=[line.split()[-1] for line in open(ProteinList,'r').readlines()]
    ListLength=len(List)
    i=0
    os.system('cp  /sfs/lustre/bahamut/scratch/jws6pq/CMfiles/SARS2.fasta  /sfs/lustre/bahamut/scratch/jws6pq/Notebook/Finished/Fastas/SARS2wS1.fasta')
    while i<ListLength:
        os.system("cat /sfs/lustre/bahamut/scratch/jws6pq/Notebook/Finished/Fastas/SARS2w" + List[i] + "S1.fasta >> /sfs/lustre/bahamut/scratch/jws6pq/Notebook/Finished/Fastas/SARS2wS1.fasta")
        i+=1
    os.system('module load gcc/9.2.0 && module load muscle/3.8.31 && muscle -in /sfs/lustre/bahamut/scratch/jws6pq/Notebook/Finished/Fastas/SARS2wS1.fasta -out /sfs/lustre/bahamut/scratch/jws6pq/Notebook/Alignment/SARS2wS1.fasta.aln')
MultipleSequenceAlignment('List')
