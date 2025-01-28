#
#  Copyright (C) 2015 Greg Landrum
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" Contains an implementation of Physicochemical property fingerprints, as
described in:
Kearsley, S. K. et al.
"Chemical Similarity Using Physiochemical Property Descriptors."
J. Chem.Inf. Model. 36, 118-127 (1996)

The fingerprints can be accessed through the following functions:
- GetBPFingerprint
- GetBTFingerprint

"""
import os.path
import re

from rdkit import Chem, RDConfig
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdFingerprintGenerator

numPathBits = rdMolDescriptors.AtomPairsParameters.numPathBits
_maxPathLen = (1 << numPathBits) - 1  # Unused variable
numFpBits = numPathBits + 2 * rdMolDescriptors.AtomPairsParameters.codeSize
fpLen = 1 << numFpBits


def _readPattyDefs(fname=os.path.join(RDConfig.RDDataDir, 'SmartsLib', 'patty_rules.txt')):
  with open(fname, 'r') as inf:
    lines = [x.strip().split('# ')[0].strip() for x in inf]
  splitl = [re.split('[ ]+', x) for x in lines if x != '']
  matchers = []
  for tpl in splitl:
    if len(tpl) > 1:
      mol = Chem.MolFromSmarts(tpl[0])
      if mol is None:
        continue
      matchers.append((mol, tpl[1]))
  return matchers


_pattyDefs = None


def AssignPattyTypes(mol, defns=None):
  """

    >>> from rdkit import Chem
    >>> AssignPattyTypes(Chem.MolFromSmiles('OCC(=O)O'))
    ['POL', 'HYD', 'OTH', 'ANI', 'ANI']

    """
  global _pattyDefs
  if defns is None:
    if _pattyDefs is None:
      _pattyDefs = _readPattyDefs()
    defns = _pattyDefs
  res = [''] * mol.GetNumAtoms()
  for matcher, nm in defns:
    matches = mol.GetSubstructMatches(matcher, uniquify=False)
    for match in matches:
      res[match[0]] = nm
  return res


typMap = dict(CAT=1, ANI=2, POL=3, DON=4, ACC=5, HYD=6, OTH=7)

_apFPG = None
_ttFPG = None


def _atomPairFingerprintFunc(mol, atomInvariants):
  global _apFPG
  if _apFPG is None:
    _apFPG = rdFingerprintGenerator.GetAtomPairGenerator()
  return _apFPG.GetSparseCountFingerprint(mol, customAtomInvariants=atomInvariants)


def _topologicalTorsionsFingerprintFunc(mol, atomInvariants):
  global _ttFPG
  if _ttFPG is None:
    _ttFPG = rdFingerprintGenerator.GetTopologicalTorsionGenerator()
  return _ttFPG.GetSparseCountFingerprint(mol, customAtomInvariants=atomInvariants)


def GetBPFingerprint(mol, fpfn=_atomPairFingerprintFunc):
  """
    >>> from rdkit import Chem
    >>> fp = GetBPFingerprint(Chem.MolFromSmiles('OCC(=O)O'))
    >>> fp.GetTotalVal()
    10
    >>> nze = fp.GetNonzeroElements()
    >>> sorted([(k, v) for k, v in nze.items()])
    [(32834, 1), (49219, 2), (98370, 2), (98401, 1), (114753, 2), (114786, 1), (114881, 1)]

    """
  typs = [typMap[x] for x in AssignPattyTypes(mol)]
  return fpfn(mol, typs)


def GetBTFingerprint(mol, fpfn=_topologicalTorsionsFingerprintFunc):
  """
    >>> from rdkit import Chem
    >>> mol = Chem.MolFromSmiles('OCC(N)O')
    >>> AssignPattyTypes(mol)
    ['POL', 'HYD', 'HYD', 'CAT', 'POL']
    >>> fp = GetBTFingerprint(mol)
    >>> fp.GetTotalVal()
    2
    >>> nze = fp.GetNonzeroElements()
    >>> sorted([(k, v) for k, v in nze.items()])
    [(538446850..., 1), (538446852..., 1)]

    """
  return GetBPFingerprint(mol, fpfn=fpfn)


# ------------------------------------
#
#  doctest boilerplate
#
def _runDoctests(verbose=None):  # pragma: nocover
  import doctest
  import sys
  failed, _ = doctest.testmod(optionflags=doctest.ELLIPSIS, verbose=verbose)
  sys.exit(failed)


if __name__ == '__main__':  # pragma: nocover
  _runDoctests()
