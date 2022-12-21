import os


def RunFoldX(Protein,Domain,ComparisonScore=827.13):
    import os
    os.system('/scratch/jws6pq/FoldX/foldx_20221231 --command=Stability --clean-mode=1 --pdb=SARS2w'+Protein+Domain+'.pdb --pdb-dir=/scratch/jws6pq/Notebook/PDB --output-dir=/scratch/jws6pq/Notebook/FoldXResults/')
def FoldXStability(Protein,Domain,ComparisonScore=827.13):
    FoldxScore=float(open('/scratch/jws6pq/Notebook/FoldXResults/SARS2w'+Protein+Domain+'_0_ST.fxout','r').readlines()[0].split()[1])
    FoldXDifference=-((ComparisonScore-FoldxScore)/ComparisonScore)*100
    return FoldXDifference

def PieceWiseRMSD(pdb_file_one, boundary_tuples, CPSplice1, CPSplice2, SpliceBoundary1, SpliceBoundary2, pdb_file_two='6VSB_B.pdb'):
    from pymol import cmd
    # ComparisonProtein=CP
    cmd.load(pdb_file_two, object='CP')
    cmd.load(pdb_file_one, object='Protein')
    cmd.remove('organic')
    cmd.select('CPNonDomain', selection='CP and not resi '+str(CPSplice1)+'-'+str(CPSplice2))
    cmd.select('CPDomain', selection='CP and resi '+str(CPSplice1)+'-'+str(CPSplice2))
    # for boundaries in boundary_tuples:

    Domainindex = pdb_file_one + ' and resi ' + str(SpliceBoundary1) + '-' + str(SpliceBoundary2)
    NonDomainindex = pdb_file_one + ' and not resi ' + str(SpliceBoundary1) + '-' + str(SpliceBoundary2)
    DomainName = pdb_file_one + 'Domain'
    NonDomainName = pdb_file_one + 'NonDomain'
    cmd.select(DomainName, selection=Domainindex)
    cmd.select(NonDomainName, selection=NonDomainindex)
    DomainRMSD=cmd.align('CPDomain', DomainName)[0]
    cmd.delete('all')
    return DomainRMSD
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

def SequenceSimilarity(Protein,Domain):
    EmbossScore=open(Protein+Domain+'.emboss','r').readlines()[25].split()[-1]
    return EmbossScore.replace('(','').replace(')','').replace('%','')
def OverallConfidence(plddt_file):
    plddt= [float(score) for score in open(plddt_file, 'r').readlines()]
    average_plddt=sum(plddt)/len(plddt)
    return average_plddt


def confidence_comparison(native_plddt, chimera_plddt, chimera_boundary_tuple, native_boundary_tuple):
    native_protein_score = [float(score) for score in open(native_plddt, 'r').readlines()]
    chimera_score=[float(score) for score in open(chimera_plddt, 'r').readlines()]
    splice_length=len(chimera_score[chimera_boundary_tuple[0]:chimera_boundary_tuple[1]])
    relative_difference=0
    native_range=[ind + native_boundary_tuple[0] for ind, x in enumerate(native_protein_score[native_boundary_tuple[0]:native_boundary_tuple[1]])]
    chimera_range = [ind + chimera_boundary_tuple[0] for ind, x in enumerate(chimera_score[chimera_boundary_tuple[0]:chimera_boundary_tuple[1]])]
    for x,y in zip(chimera_range,native_range):
        relative_difference += (chimera_score[x]-native_protein_score[y])/native_protein_score[y]*100
    relative_difference=relative_difference

    return relative_difference,splice_length

def averaging_multimer_plddt(plddt_file, subunits=3):
    """This function takes a plddt and averages the scores
    for each residue position across the number of subunints specified"""
    # Using list comprehension to turn the plddt file into a list of floats
    multimer_plddt=[float(score) for score in open(plddt_file, 'r').readlines()]
    # Calculating the length a subunits to have for step size when iterating through the list later
    monomer_length=int(len(multimer_plddt) / int(subunits))
    # creating a file to input the averaged scores
    new_plddt_file = open('Avg' + plddt_file, 'w')
    # using list comprehension to step through each the residue position of each subunit and
    # collect their scores, averaged them and return them to the new list
    averaged_scores=[sum(multimer_plddt[residue_index::monomer_length])/subunits for residue_index in range(monomer_length)]
    # Looping through the new list and inputing the averaged scores into the new file that was created
    for score in averaged_scores:
        new_plddt_file.write(str(score) + '\n')
    new_plddt_file.close()
    return new_plddt_file.name


def calc_residue_dist(residue_one, residue_two) :
    """Returns the C-alpha distance between two residues"""
    # Mostly copied from https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/protein_contact_map/
    from numpy import sqrt
    # Using distance formula to calculate residue distance
    diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord
    return sqrt(sum(diff_vector * diff_vector))
def residue_dist_matrix(pdb_file,chain_identifier,make_it_binary='Yes',distance_cutoff=7):
    """Takes a distance matrix of a protein and converts it to a binary matrix that uses 1 to signify a residue
    contact and 0 for no contact or Returns a matrix of C-alpha distances between residues in a protein chain or b"""
    # Mostly copied from https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/protein_contact_map/
    from numpy import array
    from Bio.PDB import PDBParser
    RESIDUE_NUMBER_CUTOFF=6
    # Calling the PDB structure with the parser
    STRUCTURE = PDBParser().get_structure(pdb_file, pdb_file)
    # Then selecting the chain
    PROTEIN = STRUCTURE[0][chain_identifier]
    if make_it_binary=='Yes':
        # Using list comprehension to create a nested list
        # then an array(protein length by protein length) of 1 to signifying a residue contact below distance_cutoff
        # and at least 6 residue separation or 0 for no contact
        return array([[1 if (calc_residue_dist(res_one, res_two)) <= distance_cutoff and (row-col) >= RESIDUE_NUMBER_CUTOFF else 0 for row, res_one in enumerate(PROTEIN)] for col, res_two in enumerate(PROTEIN)])
    else:
        # Using list comprehension to create a list then an array of each residue distance
        # to the rest of the chain and stacking those tuples to make a matrix
        return array([[calc_residue_dist(res_one,res_one) for res_one in PROTEIN] for res_one in PROTEIN])


def get_residue_contact_pairs(pdb_filename,chain_identifier):
    """Takes a pdb structure and return a nested list where each residue index in the list has a list with the python
    indexes for the residues it's in contact with.
    Residue indexes that have no contacts will have an empty list"""
    from numpy import where
    dist_matrix = residue_dist_matrix(pdb_filename,chain_identifier)
    x_axis,y_axis=list(where(dist_matrix==1)[0]),list(where(dist_matrix==1)[1])
    list_of_contact_pairs=[[] for rows in dist_matrix]
    for x, y in zip(x_axis, y_axis):
        list_of_contact_pairs[x].append(y)
    return list_of_contact_pairs


def correct_residue_position_for_alignment(pdb_file,chain_identifier,sequence_from_alignment):
    """Takes a nested list of residue positions
    and translates them into the residue position they are found in the alignment given"""
    contact_pairs = get_residue_contact_pairs(pdb_file,chain_identifier)
    sequence_indexing = [ind for ind, x in enumerate(sequence_from_alignment) if x != '-']
    index_dictionary = {real_index:alignment_index for real_index, alignment_index in enumerate(sequence_indexing)}
    updated_contact_map = [[index_dictionary[real_index] for real_index in contact_pairs[alignment_index]] for alignment_index, pairs in enumerate(contact_pairs)]
    return updated_contact_map
def ContactOverlap(Alignmentfile,comparison,reference='6vsb_B'):
    # Alignment in FASTA format. Make sure your benchmark sequence is first
    Sequences = open(Alignmentfile, "r").read().split('>')
    SequenceDictionary={sequence.split('\n')[0]:sequence.split('\n')[1].strip() for sequence in Sequences if len(sequence)!=0}
    ReferenceSequence,ComparisonSequence=SequenceDictionary[reference],SequenceDictionary[comparison]
    #CP is Comparison Protein and RP is Reference Protein
    CPUpdatedContactMap=CorrectResiduePositionforAlignment(comparison,ComparisonSequence)
    RPUpdatedContactMap = CorrectResiduePositionforAlignment(reference, ReferenceSequence)
    ReferenceContactMap,ComparisonContactMap=[],[]
    j=0
    #Should this be a function?
    RPContactCount,CPContactCount=0,0
    for x in ReferenceSequence:
        if x.isalpha():
            ReferenceContactMap.append([x]+RPUpdatedContactMap[j])
            RPContactCount+=len(RPUpdatedContactMap[j])+1
            j+=1
        else:
            ReferenceContactMap.append([x])
    j=0
    for x in ComparisonSequence:
        if x.isalpha():
            ComparisonContactMap.append([x]+CPUpdatedContactMap[j])
            CPContactCount+=len(CPUpdatedContactMap[j])+1
            j+=1
        else:
            ComparisonContactMap.append([x])
    TotalContacts = RPContactCount+CPContactCount
    for x,y in zip(ReferenceContactMap,ComparisonContactMap):
        for w in x:
            if w not in y and w!='-' and len(x)>1 or w not in y and w!='-' and len(y)>1:
                TotalContacts+=-1
        for v in y:
            if v not in x and v!='-' and len(x)>1 or v not in x and v!='-' and len(y)>1:
                TotalContacts += -1
    # system('/scratch/jws6pq/EMBOSS-6.6.0/emboss/needle -sprotein -gapopen 10 -gapextend 0.5 -outfile /gpfs/gpfs0/scratch/jws6pq/Notebook/Emboss/Full' + key + '.emboss -asequence /gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/' + key + '.fasta -bsequence /gpfs/gpfs0/scratch/jws6pq/BridNotebook/Fastas/SARS2.fasta')
    EmbossScore = open('/gpfs/gpfs0/scratch/jws6pq/Notebook/Emboss/Full' + comparison + '.emboss', 'r').readlines()[25].split()[-1]
    return comparison,TotalContacts,EmbossScore.replace('(','').replace(')','').replace('%','')
#Do i consider all the times where there are residues beyond SARS???????
def FaultScan(proteinpdb):
    return 1 if pymol_rmsd(proteinpdb) > 35.5 else 0
def RMSF(Protein,Timestepinps):
    from os import system as sys
    sys('echo 1 | gmx_mpi trjconv -f /gpfs/gpfs0/scratch/jws6pq/Gromacs/'+Protein+'.xtc -s /gpfs/gpfs0/scratch/jws6pq/Gromacs/3merSARS2wBatHKU4S1_production_1.tpr -dt '+Timestepinps+' -o /scratch/jws6pq/Gromacs/'+Protein+'RMSF.xtc')
    sys('echo 1 | gmx_mpi rmsf -dt '+Timestepinps+' -res -f /gpfs/gpfs0/scratch/jws6pq/Gromacs/'+Protein+'RMSF.xtc -s /gpfs/gpfs0/scratch/jws6pq/Gromacs/'+Protein+'_production_1.tpr -o '+Protein+'rmsf.xvg')