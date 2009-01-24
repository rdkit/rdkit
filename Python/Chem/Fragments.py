# $Id$
#
# Copyright (C) 2002-2006 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
""" functions to match a bunch of fragment descriptors from a file

No user-servicable parts inside.  ;-)

"""
import os
from pyRDKit import RDConfig
from pyRDKit import Chem


defaultPatternFileName = os.path.join(RDConfig.RDDataDir,'FragmentDescriptors.csv')

def _CountMatches(mol,patt,unique=1):
  res = 0
  res = len(mol.GetSubstructMatches(patt))
  return res

fns = []
def _LoadPatterns(fileName=None):
  if fileName is None:
    fileName = defaultPatternFileName
  try:
    inF = open(fileName,'r')
  except IOError:
    pass
  else:
    for line in inF.readlines():
      if len(line) and line[0] != '#':
        splitL = line.split('\t')
        if len(splitL)>=3:
          name = splitL[0]
          descr = splitL[1]
          sma = splitL[2]
          descr=descr.replace('"','')
          ok=1
          try:
            patt = Chem.MolFromSmarts(sma)
          except:
            ok=0
          else:
            if not patt or patt.GetNumAtoms()==0: ok=0
          if not ok: raise ImportError,'Smarts %s could not be parsed'%(repr(sma))
          fn = lambda x,y=1,z=patt:_CountMatches(x,z,unique=y)
          fn.__doc__ = descr
          name = name.replace('=','_')
          name = name.replace('-','_')
          fns.append((name,fn))
        
_LoadPatterns()            
for name,fn in fns:
  exec('%s=fn'%(name))
fn=None

