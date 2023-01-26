def generate_alphafold_plddt(alphafold_folder, new_plddt_file):
    """Creates a text file containing the plddt values of the highest ranked model, extracted from alphafold's result pkl file."""
    from pickle import load as p_load
    from json import load as j_load
    from os.path import exists
    from numpy import savetxt
    # Checking to see if ranking_debug.json exists. This file is the last to be output by alphafold and is a check that
    # the pkl file you want to extract from exists, as well as to avoid errors
    if exists(alphafold_folder+'ranking_debug.json'):
        # ranking_debug is also useful for determining which result pkl file is highest ranked. The model result pkls are
        # numbered by the order they are created and not their overall confidence score. The information about their rank by
        # score is found in ranking_debug
        with open(alphafold_folder+'ranking_debug.json', 'r') as json_file:
            highest_rank_model=j_load(json_file)['order'][0][0:7]
            with open(f'{alphafold_folder}result_{highest_rank_model}_multimer_v2_pred_0.pkl', 'rb') as pkl:
                prediction_data=p_load(pkl)
                savetxt(new_plddt_file, prediction_data['plddt'], fmt="%s", delimiter=" ")


def get_sequence_similarity(emboss_file):
    """Returns sequence similarity from an emboss needle file."""
    emboss_score=open(emboss_file,'r').readlines()[25].split()[-1]
    return emboss_score.replace('(','').replace(')','').replace('%','')


def overall_confidence(plddt_file):
    """Returns the average confidence score from a protein's plddt file."""
    plddt= [float(score) for score in open(plddt_file, 'r').readlines()]
    average_plddt=sum(plddt)/len(plddt)
    return average_plddt


def relative_stability(native_plddt, chimera_plddt, chimera_boundary_tuple, native_boundary_tuple):
    """Returns the relative percent difference between the two equally sized sections of plddt scores that are outlined with
    native_boundary_tuple and chimera_boundary_tuple. relative_difference=(compared value-reference value)/reference value * 100
    Native scores are assumed to be the reference value in this formula for relative difference"""
    # Pulling the plddt values as floats that start at native_boundary_tuple[0] and chimera_boundary_tuple[0], and end at
    # native_boundary_tuple[1] and chimera_boundary_tuple[1] but dont include index [1] scores.
    native_score = [float(score) for score in open(native_plddt, 'r').readlines()][native_boundary_tuple[0]:native_boundary_tuple[1]]
    chimera_score=[float(score) for score in open(chimera_plddt, 'r').readlines()][chimera_boundary_tuple[0]:chimera_boundary_tuple[1]]
    # Recording the length of the residue scores for averaging purposes later
    splice_length=len(chimera_score)
    relative_difference=sum([(chimera-native)/native*100 for native,chimera in zip(native_score,chimera_score)])
    return relative_difference,splice_length


def averaging_multimer_plddt(plddt_file, new_plddt_file,subunits):
    """This function takes a plddt and averages the scores
    for each residue position across the number of subunints specified"""
    # Using list comprehension to turn the plddt file into a list of floats
    multimer_plddt=[float(score) for score in open(plddt_file, 'r').readlines()]
    # Calculating the length a subunits to have for step size when iterating through the list later
    monomer_length=int(len(multimer_plddt) / int(subunits))
    # creating a file to input the averaged scores
    new_plddt = open(new_plddt_file, 'w')
    # using list comprehension to step through each the residue position of each subunit and
    # collect their scores, average them and return them to the new list
    averaged_scores=[sum(multimer_plddt[residue_index::monomer_length])/subunits for residue_index in range(monomer_length)]
    # Looping through the new list and inputting the averaged scores into the new file that was created
    for score in averaged_scores:
        new_plddt.write(f'{score}\n')
    new_plddt.close()
