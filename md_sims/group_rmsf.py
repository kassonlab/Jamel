#!python
#
# Simple RMSF and RMSF analysis of simulations
# Code 2022 by Peter Kasson

import glob
import json
import numpy

gmxbin = 'srun gmx'

def analyze_one(basename):
  """Runs analysis for one trajectory.
  Args: basename: basename for trajectory aligned by chain
  Rets: RMSF aligned by chain
  """
  print(basename)
  rmsf = []
  for chain_idx in range(3):
    tmpdat = numpy.loadtxt('%s_rmsf_%d.xvg' % (basename, chain_idx),
                           comments=['#', '@'])
    rmsf.append(tmpdat[:, 1]**2)
  # re-calculate RMSF
  return numpy.sqrt(numpy.mean(numpy.vstack(rmsf), axis=0))

def analyze_jointalign(basename):
  """Runs analysis for one trajectory aligned as one unit.
  Args: basename: basename for trajectory
  Rets: RMSF aligned as whole unit
  """
  tmpdat = numpy.loadtxt('%s_rmsf.xvg' % basename,
                         comments=['#', '@'])
  return tmpdat[:, 1]

if __name__ == '__main__':
  # Parse analysis of all prod.xtc files and corresponding xvg.
  # analyze aligned by unit
  rmsf_dict = {}
  for xtcfile in glob.glob('*prod.xtc'):
    keyname = xtcfile[:-4]
    rmsf_dict[keyname] = analyze_one(keyname).tolist()
  if rmsf_dict:
    with open('rmsf_vals.json', 'w') as outfile:
      json.dump(rmsf_dict, outfile)

  # analyze aligned by chain
  rmsf_dict = {}
  for xtcfile in glob.glob('*prod.xtc'):
    keyname = xtcfile[:-4]
    rmsf_dict[keyname] = analyze_jointalign(keyname).tolist()
  if rmsf_dict:
    with open('rmsf_jointvals.json', 'w') as outfile:
      json.dump(rmsf_dict, outfile)
