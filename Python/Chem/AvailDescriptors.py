## Automatically adapted for numpy.oldnumeric Jun 27, 2008 by -c

# $Id$
#
# Copyright (C) 2001-2006 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
""" constructs the list of available descriptors

"""
import types
from pyRDKit.Chem import GraphDescriptors,MolSurf,Lipinski,Fragments,Crippen,Descriptors
from pyRDKit.Chem.EState import EState_VSA
mods = [GraphDescriptors,MolSurf,EState_VSA,Lipinski,Descriptors,Crippen,Fragments]

import numpy.oldnumeric as Numeric
from pyRDKit import Chem
otherMods = [Numeric,Chem]

others = []
for mod in otherMods:
  tmp = dir(mod)
  for name in tmp:
    if name[0] != '_':
      thing = getattr(mod,name)
      if hasattr(thing,'__call__'):
        others.append(name)

descList = []
descDict = {}
for mod in mods:
  tmp = dir(mod)

  for name in tmp:
    if name[0] != '_' and name[-1] != '_' and name not in others:
      # filter out python reference implementations:
      if name[:2]=='py' and name[2:] in tmp:
        continue
      thing = getattr(mod,name)
      if hasattr(thing,'__call__'):
        descList.append((name,thing))
        descDict[name] = thing

def Desensitize():
  """ puts in all upper and all lower case versions of each
  descriptor name
  """
  ks = descDict.keys()
  for k in ks:
    fn = descDict[k]
    descDict[k.upper()]=fn
    descDict[k.lower()]=fn
  

if __name__ == '__main__':
  from pyRDKit import Chem

  m = Chem.MolFromSmiles('CCOC')
  for name,fn in descs:
    print name,fn(m)
