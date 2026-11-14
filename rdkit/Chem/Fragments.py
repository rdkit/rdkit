# Copyright (C) 2002-2010 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" functions to match a bunch of fragment descriptors from a file

The default fragment definitions ($RDBASE/Data/FragmentDescriptors.csv) are
compiled into the C++ library, and the fr_* functions here call
rdMolDescriptors.CalcFragmentDescriptor.

No user-serviceable parts inside.  ;-)

"""
import os

from rdkit import Chem, RDConfig
from rdkit.Chem import rdMolDescriptors

defaultPatternFileName = os.path.join(RDConfig.RDDataDir, 'FragmentDescriptors.csv')


def _CountMatches(mol, patt, unique=True):
  return len(mol.GetSubstructMatches(patt, uniquify=unique))


fns = []


def _makeCxxFn(name, descr):

  def fn(mol, countUnique=True):
    return rdMolDescriptors.CalcFragmentDescriptor(mol, name, countUnique=bool(countUnique))

  fn.__name__ = fn.__qualname__ = name
  fn.__doc__ = descr
  return fn


def _LoadCxxPatterns():
  for name, descr in zip(rdMolDescriptors.GetFragmentDescriptorNames(),
                         rdMolDescriptors.GetFragmentDescriptorDescriptions()):
    fns.append((name, _makeCxxFn(name, descr)))


def _LoadPatterns(fileName=None):
  """ adds the patterns in a fragment definition file to fns

  The default definitions are compiled into the C++ library, so they are not
  read from a file.

  """
  if fileName is None:
    _LoadCxxPatterns()
    return
  try:
    with open(fileName, 'r') as inF:
      for line in inF.readlines():
        if len(line) and line[0] != '#':
          splitL = line.split('\t')
          if len(splitL) >= 3:
            name = splitL[0]
            descr = splitL[1]
            sma = splitL[2]
            descr = descr.replace('"', '')
            patt = Chem.MolFromSmarts(sma)
            if not patt or patt.GetNumAtoms() == 0:
              raise ImportError('Smarts %s could not be parsed' % (repr(sma)))
            fn = lambda mol, countUnique=True, pattern=patt: _CountMatches(
              mol, pattern, unique=countUnique)
            fn.__doc__ = descr
            name = name.replace('=', '_')
            name = name.replace('-', '_')
            fns.append((name, fn))
  except IOError:
    pass


_LoadPatterns()
for name, fn in fns:
  exec('%s=fn' % (name))
fn = None
