#!/bin/bash
#SBATCH -p gpu          # partition
#SBATCH --gres=gpu:k80:1    # number of GPUs
#SBATCH -N 1            # number of nodes
#SBATCH -c 8            # number of cores
#SBATCH --mail-user=jws6pq@virginia.edu
#SBATCH --mail-type=END
#SBATCH -t 72:00:00     # time
#SBATCH -o /sfs/lustre/bahamut/scratch/jws6pq/CMfiles/slurm.out
#SBATCH -e /sfs/lustre/bahamut/scratch/jws6pq/CMfiles/slurm.err
#SBATCH -A kas_dev
#Run program
python /sfs/lustre/bahamut/scratch/jws6pq/CMfiles/alpha_list.py -s=/scratch/jws6pq/Notebook/Finished/Fastas/SARS2wTapirS1.fasta,/scratch/jws6pq/Notebook/Finished/Fastas/SARS2wSparrowHKU17S1.fasta,/scratch/jws6pq/Notebook/Finished/Fastas/TrimerSARS2.fasta -o=/scratch/jws6pq/Notebook/Finished
cd /scratch/jws6pq/Notebook/Finished
for d in ./*/ ; do (cd "$d" && python /sfs/lustre/bahamut/scratch/jws6pq/CMfiles/pickleopener.py && mydir="$(basename $PWD)" && mydirec=/scratch/jws6pq/Notebook/Plddt/$mydir && pdb=/scratch/jws6pq/Notebook/PDB/$mydir && mv overall.txt $mydir.overall && mv plddt.txt $mydirec.plddt && mv relaxed_model_1.pdb $mydir.pdb && cp $mydir.pdb /scratch/jws6pq/Notebook/PDB/ && cp $mydir.pdb /scratch/jws6pq/Notebook/Overall/ && cp $mydirec.plddt /scratch/jws6pq/Notebook/Overall/$mydir.plddt); done