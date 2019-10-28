#
# Copyright (C) 2019 Greg Landrum and T5 Informatics GmbH
#  All Rights Reserved
#
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#

from rdkit import Chem
from rdkit.Chem import rdChemReactions


class ScaffoldTreeParams(object):
  bondBreakers = ('[R:1]-!@[!R:2]>>[*:1]-[#0].[#0]-[*:2]', )
  includeGenericScaffolds = True
  includeGenericBondScaffolds = False
  keepOnlyFirstFragment = False
  keepIntermediates = True


def _addReactionsToParams(params):
  rxns = []
  for sma in params.bondBreakers:
    rxn = rdChemReactions.ReactionFromSmarts(sma)
    if rxn:
      rxns.append(rxn)
  params._breakers = tuple(rxns)


def getMolFragments(mol, params):
  """
    >>> ps = ScaffoldTreeParams()
    >>> m = Chem.MolFromSmiles('c1ccccc1CC1NC(=O)CCC1')
    >>> ms = getMolFragments(m,ps)
    >>> sorted(Chem.MolToSmiles(x) for x in ms)
  """
  if not hasattr(params, '_breakers'):
    _addReactionsToParams(params)
  res = []
  stack = [mol]
  while stack:
    wmol = stack.pop(0)
    wmol_modified = False
    for rxn in params._breakers:
      ps = rxn.RunReactants((wmol, ))
      for p in ps:
        wmol_modified = True  # we made at least one change to this
        _updateMolProps(p[0])
        stack.append(p[0])
        if params.keepIntermediates:
          res.append(p[0])
        if not params.keepOnlyFirstFragment:
          stack.append(p[1])
          if params.keepIntermediates:
            res.append(p[1])
    if not wmol_modified:
      res.append(wmol)
  return res


def makeScaffoldGeneric(mol, doBondsToo=False):
  """

    >>> m = Chem.MolFromSmiles('*c1ncc(C(=O)O)cc1')
    >>> gm = makeScaffoldGeneric(m)
    >>> Chem.SanitizeMol(gm)
    rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_NONE
    >>> Chem.MolToSmiles(gm)
    '**(=*)*1:*:*:*(*):*:*:1'
    >>> gm2 = makeScaffoldGeneric(m,doBondsToo=True)
    >>> Chem.SanitizeMol(gm2)
    rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_NONE
    >>> Chem.MolToSmiles(gm2)
    '**1***(*(*)*)**1'

    The original molecule is not affected:
    >>> Chem.MolToSmiles(m)
    '*c1ccc(C(=O)O)cn1'
    
  """
  res = Chem.Mol(mol)
  for at in res.GetAtoms():
    at.SetAtomicNum(0)
  if doBondsToo:
    for bond in res.GetBonds():
      bond.SetBondType(Chem.BondType.SINGLE)
  return res


# ------------------------------------
#
#  doctest boilerplate
#
def _runDoctests(verbose=None):  # pragma: nocover
  import sys
  import doctest
  failed, _ = doctest.testmod(optionflags=doctest.ELLIPSIS, verbose=verbose)
  sys.exit(failed)


if __name__ == '__main__':  # pragma: nocover
  _runDoctests()
