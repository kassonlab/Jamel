#!python
#
# Script to take AlphaFold output PDBs and set them up for production simulations.
# Code 2022 by Peter Kasson

import glob
import os

gmxbin = '$HOME/gromacs/build-2022/bin/gmx'
pdb2gmx = 'pdb2gmx -ff charmm36-jul2022 -water tip3p -ignh'

def prep(pdb):
  shortname = pdb[:-4]
  os.system('%s %s -f %s -o %s -p %s -i %s'
            % (gmxbin, pdb2gmx, shortname, shortname, shortname, shortname))
  os.system('%s editconf -f %s.gro -o %s_box -bt o -c -d 1'
            % (gmxbin, shortname, shortname))
  os.system('%s grompp -p %s -c %s_box -f minim.mdp -o %s_em -maxwarn 1'
            % (gmxbin, shortname, shortname, shortname))
  os.system('%s mdrun -v -deffnm %s_em' % (gmxbin, shortname))
  os.system('%s grompp -p %s -c %s -f minim.mdp -o %s_em -maxwarn 1'
            % (gmxbin, shortname, shortname, shortname))
  os.system('%s solvate -cp %s_em -o %s_sol -p %s'
            % (gmxbin, shortname, shortname, shortname))
  os.system('%s grompp -p %s -c %s_sol -f minim.mdp -o %s_preion -maxwarn 1'
            % (gmxbin, shortname, shortname, shortname))
  os.system('echo 13 | %s genion -s %s_preion -o %s_ions -p %s -neutral -conc 0.15'
            % (gmxbin, shortname, shortname, shortname))
  os.system('%s grompp -p %s -c %s_ions -f minim.mdp -o %s_em -maxwarn 1'
            % (gmxbin, shortname, shortname, shortname))
  os.system('%s mdrun -v -deffnm %s_em' % (gmxbin, shortname))
  os.system('%s grompp -p %s -c %s_em -f startup.mdp -o %s_start -maxwarn 1'
            % (gmxbin, shortname, shortname, shortname))
  os.system('%s mdrun -v -deffnm %s_start' % (gmxbin, shortname))
  os.system('%s grompp -p %s -c %s_start -f prod.mdp -o %s_prod -maxwarn 1'
            % (gmxbin, shortname, shortname, shortname))

if __name__ == '__main__':
  for pdbfile in glob.glob('3mer*pdb'):
      prep(pdbfile)
