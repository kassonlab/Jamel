def run_foldx(Protein, Domain, ComparisonScore=827.13):
    import os
    os.system('/scratch/jws6pq/FoldX/foldx_20221231 --command=Stability --clean-mode=1 --pdb=SARS2w'+Protein+Domain+'.pdb --pdb-dir=/scratch/jws6pq/Notebook/PDB --output-dir=/scratch/jws6pq/Notebook/FoldXResults/')
def get_foldx_score(Protein, Domain, ComparisonScore=827.13):
    FoldxScore=float(open('/scratch/jws6pq/Notebook/FoldXResults/SARS2w'+Protein+Domain+'_0_ST.fxout','r').readlines()[0].split()[1])
    FoldXDifference=-((ComparisonScore-FoldxScore)/ComparisonScore)*100
    return FoldXDifference


def pymol_rmsd(pdb_file_one, pdb_file_two='3merSARS2.pdb'):
    """This function takes two proteins and uses pymol to calculate the RMSD between the entirety of both proteins.
    This must be run in pymol"""
    from pymol import cmd
    cmd.load(pdb_file_two, object='CP')
    cmd.load(pdb_file_one, object='Protein')
    cmd.remove('organic')
    rmsd = cmd.align('CP', 'Protein')[0]
    cmd.delete('all')
    return rmsd


def get_sequence_similarity(emboss_file):
    emboss_score=open(emboss_file,'r').readlines()[25].split()[-1]
    return emboss_score.replace('(','').replace(')','').replace('%','')


def overall_confidence(plddt_file):
    plddt= [float(score) for score in open(plddt_file, 'r').readlines()]
    average_plddt=sum(plddt)/len(plddt)
    return average_plddt


def relative_stability(native_plddt, chimera_plddt, chimera_boundary_tuple, native_boundary_tuple):
    native_protein_score = [float(score) for score in open(native_plddt, 'r').readlines()][native_boundary_tuple[0]:native_boundary_tuple[1]]
    chimera_score=[float(score) for score in open(chimera_plddt, 'r').readlines()][chimera_boundary_tuple[0]:chimera_boundary_tuple[1]]
    splice_length=len(chimera_score)
    print(splice_length)
    relative_difference=0
    for i in range(splice_length):
        relative_difference += (chimera_score[i]-native_protein_score[i])/native_protein_score[i]*100
    return relative_difference,splice_length


relative_stability('AvgHA171to371.plddt','AvgHA171to371Spike541to741.plddt',[0, 200], [170, 370])

def averaging_multimer_plddt(plddt_file, subunits=3):
    """This function takes a plddt and averages the scores
    for each residue position across the number of subunints specified"""
    # Using list comprehension to turn the plddt file into a list of floats
    multimer_plddt=[float(score) for score in open(plddt_file, 'r').readlines()]
    # Calculating the length a subunits to have for step size when iterating through the list later
    monomer_length=int(len(multimer_plddt) / int(subunits))
    # creating a file to input the averaged scores
    new_plddt_file = open(f'Avg{plddt_file}', 'w')
    # using list comprehension to step through each the residue position of each subunit and
    # collect their scores, averaged them and return them to the new list
    averaged_scores=[sum(multimer_plddt[residue_index::monomer_length])/subunits for residue_index in range(monomer_length)]
    # Looping through the new list and inputing the averaged scores into the new file that was created
    for score in averaged_scores:
        new_plddt_file.write(f'{score}\n')
    new_plddt_file.close()
    return new_plddt_file.name


def fault_scan(proteinpdb):
    return 1 if pymol_rmsd(proteinpdb) > 35.5 else 0
