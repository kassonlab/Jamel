def call_gromacs_in_rivanna():    from os import system as syst    syst('module load gcc/9.2.0 && module load cuda/11.0.228 && module load openmpi/3.1.6 && module load '         'gcccuda/9.2.0_11.0.228 && module load goolfc/9.2.0_3.1.6_11.0.228 && module load gromacs/2021.2 && cd '         '/scratch/jws6pq/Gromacs/')def rmsd(pdb_file_one, pdb_file_two='3merSARS2.pdb'):    """This function takes two proteins and uses pymol to calculate the RMSD between the entirety of both proteins.    This must be run in pymol"""    from os import system as syst    return rmsddef truncate_xtc_file(tpr_file, xtc_file, new_xtc_file, timestep_in_ps):    from os import system as syst    syst(f'echo 1 0 | gmx_mpi trjconv -pbc cluster -f {xtc_file} -s {tpr_file} -dt {timestep_in_ps} -o {new_xtc_file}')def create_pdb_from_trajectory(tpr_file, xtc_file, new_pdb_file, timestamp_in_ps):    from os import system as syst    syst(f'echo 1 1 | gmx_mpi trjconv -f {xtc_file} -center 1 -s {tpr_file} -dump {timestamp_in_ps} -o {new_pdb_file}')    syst(f'echo 1 1 | gmx_mpi trjconv -f {new_pdb_file} -pbc mol -s {tpr_file}  -o {new_pdb_file}')def create_rmsf_file(timestep_in_ps, tpr_file, xtc_file, new_xvg_file):    from os import system as syst    syst(f'echo 3 1 | gmx_mpi rmsf -dt {timestep_in_ps} -f {xtc_file} -s {tpr_file} -o {new_xvg_file} -res -fit yes')def averaging_multimer_rmsf(xvg_file, new_xvg_file):    from itertools import groupby    rmsf_values = [tuple(value.split()) for value in open(xvg_file, 'r').readlines() if value.split()[0].isdigit()]    rmsf_values.sort(key=lambda x: int(x[0]))    sorted_values = [list(group) for key, group in groupby(rmsf_values, lambda x: x[0])]    averaged_values = [sum([float(x[-1]) for x in y]) / len(y) for y in sorted_values]    new_xvg_file = open(new_xvg_file, 'w')    for line in averaged_values:        new_xvg_file.write(str(f'{line} \n'))    new_xvg_file.close()def create_trajectory_movie_pdb(tpr_file, xtc_file, new_pdb_file):    from os import system as syst    syst(f'echo 3 1 | gmx_mpi trjconv -s {tpr_file} -f {xtc_file} -o {new_pdb_file} -fit rot+trans')