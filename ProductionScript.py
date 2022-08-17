ProteinList=open('List','r').readlines()
import numpy as np
import RBDFinder
import AccessiontoFasta
import ChimeraGenerator
import os
NumberofProteinsPerSlurm=6
ChimeraOnly='No'
FastasforRun=np.empty(len(ProteinList)*2, dtype=object)
k=0
for line in ProteinList:
    ProteinInfo=line.split()
    # AccessiontoFasta.AccessionNumbertoFasta(ProteinInfo[-1],ProteinInfo[0])
    # AccessiontoFasta.FastatoAlignmentFinder(ProteinInfo[-1])
    SpliceBoundary=RBDFinder.AlignmentFinder('/sfs/lustre/bahamut/scratch/jws6pq/Notebook/Alignment/'+ProteinInfo[-1]+'onSARS2.aln',SequenceofInterest='AYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLSEFRVYSSANNCTFEYVSQPFLKNLREFVFKNIDGYFKIYSKHTPPQGFSALEPLVDLPIGINITRFQTLLAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSYNYLYRNLKPFERDISTEIYNCYFPLQSYGFQPTVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICA')
    FastasforRun[k]=ChimeraGenerator.DomainExchange('SARS2.fasta',ProteinInfo[-1]+'.fasta',SpliceBoundary[0],SpliceBoundary[1],Boundary1=1,Boundary2=540,Domain='S1')[0]
    FastasforRun[k+1]=ChimeraGenerator.DomainExchange('SARS2.fasta',ProteinInfo[-1]+'.fasta',SpliceBoundary[0],SpliceBoundary[1],Boundary1=1,Boundary2=540,Domain='S1')[1]
    k+=2
i=0
j=0
k=0
while i in range(len(FastasforRun)-NumberofProteinsPerSlurm-1):
    os.system('cp /sfs/lustre/bahamut/scratch/jws6pq/CMfiles/MultimerAlphaFold.slurm /sfs/lustre/bahamut/scratch/jws6pq/CMfiles/'+str(j)+'MultimerAlphaFold.slurm')
    file=open('/sfs/lustre/bahamut/scratch/jws6pq/CMfiles/'+str(j)+'MultimerAlphaFold.slurm','a')
    if ChimeraOnly=='No':
        file.write('python /sfs/lustre/bahamut/scratch/jws6pq/CMfiles/multimeralpha_list.py -s='+FastasforRun[i]+FastasforRun[i+1]+FastasforRun[i+2]+FastasforRun[i+3]+FastasforRun[i+4]+FastasforRun[i+5]+' -o=/scratch/jws6pq/Notebook/Finished\n')
        i+=NumberofProteinsPerSlurm
    elif ChimeraOnly=='Yes':
        file.write('python /sfs/lustre/bahamut/scratch/jws6pq/CMfiles/multimeralpha_list.py -s='+FastasforRun[i+1]+FastasforRun[i+3]+FastasforRun[i+5]+FastasforRun[i+7]+FastasforRun[i+9]+FastasforRun[i+11]+' -o=/scratch/jws6pq/Notebook/Finished\n')
        i+=NumberofProteinsPerSlurm*2
    # file.write('cd /scratch/jws6pq/Notebook/Finished\n')
    # file.write('for d in ./*/ ; do (cd "$d" && python /sfs/lustre/bahamut/scratch/jws6pq/CMfiles/pickleopener.py && mydir="$(basename $PWD)" && mydirec=/scratch/jws6pq/Notebook/Plddt/$mydir && pdb=/scratch/jws6pq/Notebook/PDB/$mydir && mv overall.txt $mydir.overall && mv plddt.txt $mydirec.plddt && mv ranked_0.pdb $mydir.pdb && cp $mydir.pdb /scratch/jws6pq/Notebook/PDB/ && cp $mydir.pdb /scratch/jws6pq/Notebook/Overall/ && cp $mydirec.plddt /scratch/jws6pq/Notebook/Overall/$mydir.plddt); done')
    file.close()
    # os.system('sbatch /sfs/lustre/bahamut/scratch/jws6pq/CMfiles/'+str(j)+'MultimerAlphaFold.slurm')
    j+=1