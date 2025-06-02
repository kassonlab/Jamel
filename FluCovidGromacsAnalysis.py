from Chimeragenesis.Scripts import GromacsAnalysis
from json import load
from sys import argv
from numpy import savetxt,array
argument_json=str(argv[1])
protein_list=str(argv[2])
subunits=int(argv[3])

with open(argument_json, 'rb') as jfile:
    argument_dict=load(jfile)["arguments"]
with open(protein_list, 'r') as p_list:
    p_list=p_list.readlines()
    protein_list= [x.split()[0] for x in p_list]

start_time=argument_dict['rmsd_timestep_ps'][0]
stop_time=argument_dict['rmsd_timestep_ps'][1]+1
rmsd_step=argument_dict['rmsd_timestep_ps'][2]
character_to_replace =argument_dict["character_to_replace"]
tpr_file=[argument_dict['tpr_naming'].replace(character_to_replace, protein) for protein in protein_list]
xtc_file=[argument_dict['xtc_naming'].replace(character_to_replace, protein) for protein in protein_list]
time_step=[argument_dict['rmsf_timestep_ps'] for protein in protein_list]
new_xtc_file=[argument_dict['truncated_xtc_naming'].replace(character_to_replace,protein) for protein in protein_list]
xvg_file=[argument_dict['rmsf_naming'].replace(character_to_replace,protein) for protein in protein_list]
movie_pdb=[argument_dict['movie_pdb_file'].replace(character_to_replace,protein) for protein in protein_list]
rmsd_time_steps=range(start_time,stop_time,rmsd_step)

if argument_dict['truncated_xtc_naming'][0]!='#':
    list(map(GromacsAnalysis.truncate_xtc_file, tpr_file, xtc_file, new_xtc_file, time_step))
if argument_dict['rmsf_naming'][0] != '#':
    list(map(GromacsAnalysis.create_rmsf_file, time_step, tpr_file, new_xtc_file, xvg_file))
if argument_dict['movie_pdb_file'][0] != '#':
    list(map(GromacsAnalysis.create_trajectory_movie_pdb, tpr_file, new_xtc_file, movie_pdb))
if subunits>1 and argument_dict['averaged_rmsf_naming'][0]!='#':
    averaged_xvg_file=[argument_dict['averaged_rmsf_naming'].replace(character_to_replace,protein) for protein in protein_list]
    list(map(GromacsAnalysis.averaging_multimer_rmsf, xvg_file, averaged_xvg_file))
if argument_dict['rmsd_naming'][0] != '#':
    rmsd_change=[]
    for protein_index, protein in enumerate(protein_list):
        rmsd_files=[argument_dict['rmsd_naming'].replace(character_to_replace,protein,1).replace(character_to_replace,str(time)) for time in rmsd_time_steps]
        [GromacsAnalysis.create_pdb_from_trajectory(tpr_file[protein_index], new_xtc_file[protein_index], file, time) for file,time in zip(rmsd_files, rmsd_time_steps)]
        rmsd=[GromacsAnalysis.pymol_rmsd(rmsd_files[0], file) for file in rmsd_files]
        rmsd_change.append([protein]+rmsd)
    savetxt(argument_dict['rmsd_results_csv'], array(rmsd_change, dtype=object),
            fmt=','.join(['%s']+['%s' for time in rmsd_time_steps]))
