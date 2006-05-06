# $Id: Descriptors.py 5007 2006-02-22 15:14:41Z glandrum $
#
# Copyright (C) 2001-2006 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
import Chem

def MolWt(mol,heavyAtomsOnly=0):
  """ The average molecular weight of the molecule

  >>> '%.2f'%MolWt(Chem.MolFromSmiles('CC'))
  '30.07'
  >>> '%.2f'%MolWt(Chem.MolFromSmiles('CC'),1)
  '24.02'
  >>> '%.2f'%MolWt(Chem.MolFromSmiles('[NH4+].[Cl-]'))
  '53.49'
  """
  hMass = Chem.GetPeriodicTable().GetAtomicWeight(1)
  accum = 0.0
  for atom in mol.GetAtoms():
    accum += atom.GetMass()
    if not heavyAtomsOnly:
      accum += atom.GetTotalNumHs()*hMass
  return accum
MolWt.version="1.0.0"

HeavyAtomMolWt=lambda x:MolWt(x,1)
HeavyAtomMolWt.__doc__="""The average molecular weight of the molecule ignoring hydrogens

  >>> '%.2f'%HeavyAtomMolWt(Chem.MolFromSmiles('CC'))
  '24.02'
  >>> '%.2f'%HeavyAtomMolWt(Chem.MolFromSmiles('[NH4+].[Cl-]'))
  '49.46'

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
  return doctest.testmod(sys.modules["__main__"])

if __name__ == '__main__':
  import sys
  failed,tried = _test()
  sys.exit(failed)
