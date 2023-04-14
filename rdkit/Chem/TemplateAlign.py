# $Id$
#
#  Copyright (C) 2005-2006 Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
from rdkit import Chem, Geometry
from rdkit.Chem import rdDepictor


def AlignMolToTemplate2D(
  mol,
  template,
  match=None,
  clearConfs=False,
  templateConfId=-1,
):
  """
   Arguments:

     - mol:      the molecule to be aligned
     - template: the template to align to
     - match:    If provided, this should be a sequence of integers
                 containing the indices of the atoms in mol that match
                 those in template. This is the result of calling:
                   mol.GetSubstructMatch(template)
     - clearConfs: toggles removing any existing conformers on mol

    Returns the confId of the conformer containing the depiction

    >>> patt = Chem.MolFromSmiles('C1CC1')
    >>> rdDepictor.Compute2DCoords(patt)
    0
    >>> mol = Chem.MolFromSmiles('OC1CC1CC1CCC1')
    >>> rdDepictor.Compute2DCoords(mol)
    0
    >>> pc = patt.GetConformer(0)
    >>> mc = mol.GetConformer(0)

    We start out with the molecules not aligned:

    >>> vs = [abs(pc.GetAtomPosition(i).x-mc.GetAtomPosition(i+1).x) for i in range(pc.GetNumAtoms())]
    >>> [x<1e-4 for x in vs]
    [False, False, False]

    But then we can replace the conformer of mol:

    >>> AlignMolToTemplate2D(mol,patt,clearConfs=True)
    0
    >>> mol.GetNumConformers()
    1
    >>> pc = patt.GetConformer(0)
    >>> mc = mol.GetConformer(0)
    >>> vs = [abs(pc.GetAtomPosition(i).x-mc.GetAtomPosition(i+1).x) for i in range(pc.GetNumAtoms())]
    >>> [x<1e-4 for x in vs]
    [True, True, True]

    If we like, we can specify the atom map explicitly in order to align to the second
    matching ring in the probe molecule:

    >>> match = (5,6,7)
    >>> AlignMolToTemplate2D(mol,patt,clearConfs=True,match=match)
    0
    >>> mol.GetNumConformers()
    1
    >>> pc = patt.GetConformer(0)
    >>> mc = mol.GetConformer(0)
    >>> vs = [abs(pc.GetAtomPosition(i).x-mc.GetAtomPosition(i+1).x) for i in range(pc.GetNumAtoms())]
    >>> [x<1e-4 for x in vs]
    [False, False, False]
    >>> vs = [abs(pc.GetAtomPosition(i).x-mc.GetAtomPosition(i+5).x) for i in range(pc.GetNumAtoms())]
    >>> [x<1e-4 for x in vs]
    [True, True, True]



  """
  if not match:
    match = mol.GetSubstructMatch(template)
  if not match:
    raise ValueError('no match between mol and template')

  atomMap = {}
  templateConf = template.GetConformer(templateConfId)
  for i, idx in enumerate(match):
    p = templateConf.GetAtomPosition(i)
    atomMap[idx] = Geometry.Point2D(p.x, p.y)
  molConfId = rdDepictor.Compute2DCoords(mol, clearConfs=clearConfs, coordMap=atomMap)
  return molConfId


#------------------------------------
#
#  doctest boilerplate
#
def _test():
  import doctest
  import sys

  from rdkit.Chem import rdDepictor
  rdDepictor.SetPreferCoordGen(False)
  return doctest.testmod(sys.modules["__main__"])


if __name__ == '__main__':
  import sys
  failed, tried = _test()
  sys.exit(failed)
