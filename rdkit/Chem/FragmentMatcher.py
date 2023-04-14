#
# Copyright (C) 2002-2006 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""Exposes a class for matching fragments of molecules.

The class exposes a simple API:

If you want a matcher that hits C=O, you'd do:

>>> p = FragmentMatcher()
>>> p.Init('C=O')

you can then match with:

>>> mol = Chem.MolFromSmiles('CC(=O)O')
>>> p.HasMatch(mol)
1
>>> p.HasMatch(Chem.MolFromSmiles('CC(C)C'))
0

information about the matches:

>>> len(p.GetMatches(Chem.MolFromSmiles('CC=O')))
1
>>> len(p.GetMatches(Chem.MolFromSmiles('O=CC=O')))
2

or, you can add exclusion fragments (defined as smarts) with:

>>> p.AddExclusion('c1ccccc1')

now the matcher will not hit anything that has a benzene ring.

>>> p.HasMatch(Chem.MolFromSmiles('CC=O'))
1
>>> p.HasMatch(Chem.MolFromSmiles('c1ccccc1CC=O'))
0


"""
from rdkit import Chem


class FragmentMatcher(object):

  def __init__(self):
    self._onPatt = None
    self._offPatts = []

  def AddExclusion(self, sma):
    self._offPatts.append(Chem.MolFromSmarts(sma))

  def Init(self, sma):
    self._onPatt = Chem.MolFromSmarts(sma)

  def GetSMARTS(self):
    pass

  def GetExclusionSMARTS(self):
    pass

  def HasMatch(self, mol):
    if self._onPatt is None:
      return 0
    t = mol.HasSubstructMatch(self._onPatt)
    if not t:
      return 0
    else:
      for patt in self._offPatts:
        if mol.HasSubstructMatch(patt):
          return 0
    return 1

  def GetMatch(self, mol):
    if self._onPatt is None:
      return None
    return mol.GetSubstructMatch(self._onPatt)

  def GetMatches(self, mol, uniquify=1):
    if self._onPatt is None:
      return None
    return mol.GetSubstructMatches(self._onPatt, uniquify=uniquify)

  def GetBond(self, idx):
    if self._onPatt is None:
      return None
    return self._onPatt.GetBondWithIdx(idx)


  #------------------------------------
  #
  #  doctest boilerplate
  #
def _test():
  import doctest
  import sys
  return doctest.testmod(sys.modules["__main__"])


if __name__ == '__main__':
  import sys
  failed, tried = _test()
  sys.exit(failed)
