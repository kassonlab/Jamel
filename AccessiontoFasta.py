def AccessionNumbertoFasta(Protein,Accession):
    import numpy as np
    from Bio import Entrez
    handle = Entrez.efetch(db='protein', id=Accession,retmode='text',rettype='fasta').readlines()
    #record = Entrez.read(handle)
    np.savetxt('/sfs/lustre/bahamut/scratch/jws6pq/Notebook/Finished/Fastas/'+Protein+'.fasta', handle, fmt="%s", delimiter="")
    return Protein
def FastatoAlignmentFinder(Protein):
    os.system("cp  /sfs/lustre/bahamut/scratch/jws6pq/Notebook/Finished/FastasAvianD274.fasta  /sfs/lustre/bahamut/scratch/jws6pq/Notebook/Finished/Fastas/"+Protein+"onSARS2.fasta")
    os.system("cat /sfs/lustre/bahamut/scratch/jws6pq/CMfiles/SARS2.fasta >>  /sfs/lustre/bahamut/scratch/jws6pq/Notebook/Finished/Fastas/"+Protein+"onSARS2.fasta")
    os.system('/sfs/lustre/bahamut/scratch/jws6pq/clustalw-2.1/src/bin/clustalw2  /sfs/lustre/bahamut/scratch/jws6pq/Notebook/Finished/Fastas/'+Protein+"onSARS2.fasta")

