os.system('/sfs/lustre/bahamut/scratch/jws6pq/EMBOSS-6.6.0/emboss/needle -sprotein -gapopen 10 -gapextend 0.5 -outfile /sfs/lustre/bahamut/scratch/jws6pq/Notebook/Emboss/' + Protein + '.emboss -asequence asis:' + SequenceofInterest + ' -bsequence asis:' + FoundAlignment)
open(Protein,'r')


def Analysis():
    import RBDFinder
    import PlddtComparison
    import PieceRMSD
