#!python
#
# Simple RMSD and RMSF analysis of simulations
# Code 2022 by Peter Kasson

import glob
import json
import os
import numpy

gmxbin = 'srun gmx'

def analyze_one(basename):
  """Runs analysis for one trajectory.
  Args: basename: basename for trajectory
  Rets: RMSD aligned by chain
  """
  os.system('rm test.csv*ndx')
  rmsd = []
  for chain_idx in range(3):
    os.system('%s select -s %s -on test.csv%d '
              '-select "molindex %d; molindex %d and name CA"'
              % (gmxbin, basename, chain_idx, chain_idx+1, chain_idx+1))
    os.system('echo %s| %s rms -f %s.xtc -s %s -dt 1000 -n test.csv%d -o %s_rmsd_%d'
              % ('1 1', gmxbin, basename, basename, chain_idx, basename, chain_idx))
    os.system('echo %s |%s rmsf -s %s -f %s.xtc -dt 1000 -fit -n test.csv%d -o %s_rmsf_%d -res'
              % ('1 1', gmxbin, basename, basename, chain_idx, basename, chain_idx))
    tmpdat = numpy.loadtxt('%s_rmsd_%d.xvg' % (basename, chain_idx),
                           comments=['#','@'])
    if rmsd:
      rmsd[:, 1] += tmpdat[:, 1]**2
    else:
      rmsd = tmpdat**2
  return numpy.sqrt(rmsd)

if __name__ == '__main__':
  rmsd_dict = {}
  for xtcfile in glob.glob('*prod.xtc'):
    keyname = xtcfile[:-4]
    rmsd_dict[keyname] = analyze_one(keyname).tolist()
  with open('rmsd_vals.json', 'w') as outfile:
    json.dump(rmsd_dict, outfile)
