ProteinList=open('List','r').readlines()
import numpy as np
import RBDFinder
import AccessiontoFasta
import RBDExchanger
import os
FastasforRun=np.empty(len(ProteinList), dtype=object)
k=0
for line in ProteinList:
    x=line.split()
    AccessiontoFasta.AccessionNumbertoFasta(x[-1],x[0])
    AccessiontoFasta.FastatoAlignmentFinder(x[-1])
    SpliceBoundary=RBDFinder.AlignmentFinder(x[-1]+'onSARS2.aln')
    FastasforRun[k]=RBDExchanger.RBDExchange('SARS2.fasta',x[-1]+'.fasta',SpliceBoundary[0],SpliceBoundary[1])
    k+=1
i=0
j=0
for j in range(int(len(ProteinList)/4)):
    os.system('cp /sfs/lustre/bahamut/scratch/jws6pq/CMfiles/NeoAlphaFold.slurm /sfs/lustre/bahamut/scratch/jws6pq/CMfiles/'+str(j)+'AlphaFold.slurm')
    file=open('/sfs/lustre/bahamut/scratch/jws6pq/CMfiles/'+str(j)+'AlphaFold.slurm','a')
    file.write('python /sfs/lustre/bahamut/scratch/jws6pq/CMfiles/alpha_list.py -s='+FastasforRun[i]+','+FastasforRun[i+1]+','+FastasforRun[i+2]+','+FastasforRun[i+3]+' -o=/scratch/jws6pq/Notebook/Finished\n')
    file.write('cd /scratch/jws6pq/Notebook/Finished\n')
    file.write('for d in ./*/ ; do (cd "$d" && python /sfs/lustre/bahamut/scratch/jws6pq/CMfiles/pickleopener.py && mydir="$(basename $PWD)" && mydirec=/scratch/jws6pq/Notebook/Plddt/$mydir && mv overall.txt $mydir.overall && mv plddt.txt $mydirec.plddt && mv relaxed_model_1.pdb $mydir.pdb); done')
    file.close()
    os.system('sbatch /sfs/lustre/bahamut/scratch/jws6pq/CMfiles/'+str(j)+'AlphaFold.slurm')
    i+=4