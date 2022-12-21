def AccessionNumbertoFasta(protein, fasta_file_name, accession, subunits=3):
    from Bio import Entrez
    Entrez.email = 'jws6pq@virginia.edu'
    handle = Entrez.efetch(db='protein', id=accession, retmode='text', rettype='fasta').readlines()
    sequence=''.join([x for x in handle if  x[0]!='>'  if x!='']).rstrip().strip().replace('\n','').replace(' ','')
    monomer_file=open('/gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/' + protein + '.fasta', 'w')
    monomer_file.write('>' + protein + '\n' + sequence + '\n')
    multimer_file = open('/gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/' + str(subunits) + 'mer' + fasta_file_name + '.fasta', 'w')
    for x in range(subunits):
        multimer_file.write('>' + protein + '\n' + sequence + '\n')
    multimer_file.close()
def FastatoAlignmentFinder(protein):
    from os import system
    system('cp  /scratch/jws6pq/CMfiles/SARS2.fasta  /gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/' + protein + "onSARS2.fasta")
    system("cat /gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/" + protein + ".fasta >> /gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/" + protein + "onSARS2.fasta")
    system('module load gcc/9.2.0 && module load muscle/3.8.31 && muscle -in /gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/' + protein + 'onSARS2.fasta -clw -out /scratch/jws6pq/Notebook/Alignment/' + protein + 'onSARS2.aln')
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
def pdb2fasta(pdb_file):
    from Bio import SeqIO
    with open(pdb_file, 'r') as pdb:
        for record in SeqIO.parse(pdb, 'pdb-atom'):
            print('>' + record.id)
            print(record.seq)
