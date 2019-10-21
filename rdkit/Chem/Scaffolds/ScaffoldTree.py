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


class ScaffoldTreeParams(object):
  bondsToBreak = ('[R]!@[!R]', )


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
