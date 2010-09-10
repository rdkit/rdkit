# $Id$
#
# Copyright (C) 2001-2006 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def _pyMolWt(mol,heavyAtomsOnly=0):
  """ DEPRECATED 
  """
  hMass = Chem.GetPeriodicTable().GetAtomicWeight(1)
  accum = 0.0
  for atom in mol.GetAtoms():
    accum += atom.GetMass()
    if not heavyAtomsOnly:
      accum += atom.GetTotalNumHs()*hMass
  return accum
_pyMolWt.version="1.0.0"

MolWt = lambda *x,**y:rdMolDescriptors._CalcMolWt(*x,**y)
MolWt.version=rdMolDescriptors._CalcMolWt_version
MolWt.__doc__="""The average molecular weight of the molecule ignoring hydrogens

  >>> MolWt(Chem.MolFromSmiles('CC'))
  30.07...
  >>> MolWt(Chem.MolFromSmiles('[NH4+].[Cl-]'))
  53.49...

"""

HeavyAtomMolWt=lambda x:MolWt(x,1)
HeavyAtomMolWt.__doc__="""The average molecular weight of the molecule ignoring hydrogens

  >>> HeavyAtomMolWt(Chem.MolFromSmiles('CC'))
  24.02...
  >>> HeavyAtomMolWt(Chem.MolFromSmiles('[NH4+].[Cl-]'))
  49.46...

"""
HeavyAtomMolWt.version="1.0.0"

def NumValenceElectrons(mol):
  """ The number of valence electrons the molecule has

  >>> NumValenceElectrons(Chem.MolFromSmiles('CC'))
  14.0
  >>> NumValenceElectrons(Chem.MolFromSmiles('C(=O)O'))
  18.0
  >>> NumValenceElectrons(Chem.MolFromSmiles('C(=O)[O-]'))
  18.0
  >>> NumValenceElectrons(Chem.MolFromSmiles('C(=O)'))
  12.0
  

  """
  tbl = Chem.GetPeriodicTable()
  accum = 0.0
  for atom in mol.GetAtoms():
    accum += tbl.GetNOuterElecs(atom.GetAtomicNum())
    accum -= atom.GetFormalCharge()
    accum += atom.GetTotalNumHs()

  return accum
NumValenceElectrons.version="1.0.0"

#------------------------------------
#
#  doctest boilerplate
#
def _test():
  import doctest,sys
  return doctest.testmod(sys.modules["__main__"],optionflags=doctest.ELLIPSIS)

if __name__ == '__main__':
  import sys
  failed,tried = _test()
  sys.exit(failed)
