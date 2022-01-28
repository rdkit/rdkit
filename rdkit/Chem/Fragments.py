# $Id$
#
# Copyright (C) 2002-2010 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" functions to match a bunch of fragment descriptors from a file

No user-servicable parts inside.  ;-)

"""
import os
from rdkit import RDConfig
from rdkit import Chem

defaultPatternFileName = os.path.join(RDConfig.RDDataDir, 'FragmentDescriptors.csv')


def _CountMatches(mol, patt, unique=True):
  return len(mol.GetSubstructMatches(patt, uniquify=unique))


fns = []


def _LoadPatterns(fileName=None):
  if fileName is None:
    fileName = defaultPatternFileName
  try:
    with open(fileName, 'r') as inF:
      for line in inF.readlines():
        if len(line) and line[0] != '#':
          splitL = line.split('\t')
          if len(splitL) >= 3:
            name = splitL[0].replace('=', '_').replace('-', '_')
            descr = splitL[1].replace('"', '')
            sma = splitL[2]
            patt = Chem.MolFromSmarts(sma)
            if not patt or patt.GetNumAtoms() == 0:
              raise ImportError(f'Smarts {repr(sma)} could not be parsed')
            
            fn = lambda mol, countUnique=True, pattern=patt: _CountMatches(mol, pattern, unique=countUnique)
            fn.__doc__ = descr
            fns.append((name, fn))
  except IOError:
    pass


_LoadPatterns()
for name, fn in fns:
  exec(f'{name}=fn')
fn = None
