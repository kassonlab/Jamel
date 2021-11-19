ProteinList=open(List,'r').readlines()
import RBDFinder
import AccessiontoFasta
import RBDExchanger
FastasforRun=np.empty(len(FastaList))
i=0
for line in ProteinList:
    AccessiontoFasta.AccessionNumbertoFasta(ProteinList[i,1],ProteinList[i,0])
    AccessiontoFasta.FastatoAlignmentFinder(ProteinList[i,1])
    RBDFinder.AlignmentFinder(ProteinList[i,1]+'onSARS2.aln')
    SpliceBoundary=RBDFinder.AlignmentFinder
    FastasforRun[k]=RBDExchanger.RBDExchange('SARS2.fasta',Protein+'.fasta',SpliceBoundary[0],SpliceBoundary[1])
file=open('AlphaFold.slurm','a')
file.writelines('python alpha_list.py -s= -o=/scratch/jws6pq/Notebook/Finished','cd /scratch/jws6pq/Notebook/Finished','for d in ./*/ ; do (cd "$d" && python /sfs/lustre/bahamut/scratch/jws6pq/CMfiles/pickleopener.py && mydir="$(basename $PWD)" && mv overall.txt $mydir.overall && mv plddt.txt $mydir.plddt && mv relaxed_model_1.pdb $mydir.pdb); done')
print(file)
