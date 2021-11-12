def AccessionNumbertoFasta(List)
    import numpy as np
    from Bio import Entrez
    
    Entrez.email = "jws6pw@virginia.edu"
    handle = Entrez.efetch(db='protein', id='YP_005352871',retmode='text',rettype='fasta').readlines()
    #record = Entrez.read(handle)
    np.savetxt('Wigeon.fasta', handle, fmt="%s", delimiter="")

