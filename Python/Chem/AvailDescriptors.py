# $Id$
#
# Copyright (C) 2001-2006 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
""" constructs the list of available descriptors

"""
import types
from Chem import GraphDescriptors,MolSurf,Lipinski,Fragments,Crippen,Descriptors
from Chem.EState import EState_VSA
mods = [GraphDescriptors,MolSurf,EState_VSA,Lipinski,Descriptors,Crippen,Fragments]

import Numeric,Chem
otherMods = [Numeric,Chem]

others = []
for mod in otherMods:
  tmp = dir(mod)
  for name in tmp:
    if name[0] != '_':
      thing = getattr(mod,name)
      if type(thing)==types.FunctionType:
        others.append(name)

descList = []
descDict = {}
for mod in mods:
  tmp = dir(mod)

  for name in tmp:
    if name[0] != '_' and name[-1] != '_' and name not in others:
      thing = getattr(mod,name)
      if type(thing)==types.FunctionType:
        # we need to store the name too, just in case
        #  the function is a lambda
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
  import Chem

  m = Chem.MolFromSmi('CCOC')
  for name,fn in descs:
    print name,fn(m)
