def AccessionNumbertoFasta(protein, fasta_file_name, accession, subunits=3, *mulitmer_name):
    from Bio import Entrez
    Entrez.email = 'jws6pq@virginia.edu'
    handle = Entrez.efetch(db='protein', id=accession, retmode='text', rettype='fasta').readlines()
    sequence=''.join([x for x in handle if  x[0]!='>'  if x!='']).strip().replace('\n','')
    monomer_file=open(fasta_file_name, 'w')
    title=f'>{protein}\n{sequence}\n'
    monomer_file.write(title)
    monomer_file.close()
    if subunits != 1:
        multimer_file = open(mulitmer_name, 'w')
        for x in range(subunits):
            multimer_file.write(title)
        multimer_file.close()
def FastatoAlignmentFinder(protein):
    from os import system
    system('cp  /scratch/jws6pq/CMfiles/SARS2.fasta  /gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/' + protein + "onSARS2.fasta")
    system("cat /gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/{0}.fasta >> /gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/{0}onSARS2.fasta".format(protein))
    system('module load gcc/9.2.0 && module load muscle/3.8.31 && muscle -in /gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/{0}onSARS2.fasta -clw -out /scratch/jws6pq/Notebook/Alignment/{0}onSARS2.aln'.format(protein))
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
