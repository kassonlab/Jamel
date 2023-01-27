import GromacsAnalysis
from json import load
from sys import argv
argument_json=str(argv[1])
protein_list=str(argv[2])
subunits=int(argv[3])
with open(argument_json, 'rb') as jfile:
    argument_dict=load(jfile)["arguments"]
with open(protein_list, 'r') as p_list:
    protein_list= [x.split()[0] for x in p_list.readlines()]
character_to_replace =argument_dict["character_to_replace"]
tpr_file=[argument_dict['tpr_naming'].replace(character_to_replace, protein) for protein in protein_list]
xtc_file=[argument_dict['xtc_naming'].replace(character_to_replace, protein) for protein in protein_list]
time_step=[argument_dict['rmsf_timestep_ps'] for protein in protein_list]
new_xtc_file=[argument_dict['truncated_xtc_naming'].replace(character_to_replace,protein) for protein in protein_list]
xvg_file=[argument_dict['xvg_naming'].replace(character_to_replace,protein) for protein in protein_list]
movie_pdb=[argument_dict['movie_pdb_file'].replace(character_to_replace,protein) for protein in protein_list]
rmsd_time_steps=[str(time) for time in range(0,1001,100)]
GromacsAnalysis.call_gromacs_in_rivanna()
list(map(GromacsAnalysis.truncate_xtc_file,tpr_file,xtc_file,new_xtc_file,time_step))
list(map(GromacsAnalysis.create_rmsf_file,time_step,tpr_file,new_xtc_file,xvg_file))
list(map(GromacsAnalysis.create_trajectory_movie_pdb, tpr_file, new_xtc_file, movie_pdb))
if subunits>1:
    averaged_xvg_file=[argument_dict['averaged_xvg_naming'].replace(character_to_replace,protein) for protein in protein_list]
    list(map(GromacsAnalysis.averaging_multimer_rmsf, xvg_file, averaged_xvg_file))
