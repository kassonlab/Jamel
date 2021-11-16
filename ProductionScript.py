ProteinList=open(List,'r').readlines()
import RBDFinder
import AccessiontoFasta
i=0
for line in ProteinList:
    AccessiontoFasta.AccessionNumbertoFasta(ProteinList[i,1],ProteinList[i,0])
    AccessiontoFasta.FastatoAlignmentFinder(ProteinList[i,1])
    RBDFinder.AlignmentFinder(ProteinList[i,1]+'onSARS2.aln')
