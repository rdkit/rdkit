# $Id$
#
# Copyright (C) 2001-2010 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
from rdkit import Chem
from rdkit.Chem import rdPartialCharges
import collections

def _isCallable(thing):
  return (hasattr(collections,'Callable') and isinstance(thing,collections.Callable)) or \
              hasattr(thing,'__call__')   

_descList=[]
def _setupDescriptors(namespace):
  global _descList,descList
  from rdkit.Chem import GraphDescriptors,MolSurf,Lipinski,Fragments,Crippen
  from rdkit.Chem.EState import EState_VSA
  mods = [GraphDescriptors,MolSurf,EState_VSA,Lipinski,Crippen,Fragments]

  otherMods = [Chem]

  for nm,thing in namespace.items():
    if nm[0]!='_' and _isCallable(thing):
      _descList.append((nm,thing))
  
  others = []
  for mod in otherMods:
    tmp = dir(mod)
    for name in tmp:
      if name[0] != '_':
        thing = getattr(mod,name)
        if _isCallable(thing):
          others.append(name)

  for mod in mods:
    tmp = dir(mod)

    for name in tmp:
      if name[0] != '_' and name[-1] != '_' and name not in others:
        # filter out python reference implementations:
        if name[:2]=='py' and name[2:] in tmp:
          continue
        thing = getattr(mod,name)
        if _isCallable(thing):
          namespace[name]=thing
          _descList.append((name,thing))
  descList=_descList

from rdkit.Chem import rdMolDescriptors as _rdMolDescriptors
MolWt = lambda *x,**y:_rdMolDescriptors._CalcMolWt(*x,**y)
MolWt.version=_rdMolDescriptors._CalcMolWt_version
MolWt.__doc__="""The average molecular weight of the molecule

  >>> MolWt(Chem.MolFromSmiles('CC'))
  30.07
  >>> MolWt(Chem.MolFromSmiles('[NH4+].[Cl-]'))
  53.49...

"""

HeavyAtomMolWt=lambda x:MolWt(x,True)
HeavyAtomMolWt.__doc__="""The average molecular weight of the molecule ignoring hydrogens

  >>> HeavyAtomMolWt(Chem.MolFromSmiles('CC'))
  24.02...
  >>> HeavyAtomMolWt(Chem.MolFromSmiles('[NH4+].[Cl-]'))
  49.46

"""
HeavyAtomMolWt.version="1.0.0"

ExactMolWt = lambda *x,**y:_rdMolDescriptors.CalcExactMolWt(*x,**y)
ExactMolWt.version=_rdMolDescriptors._CalcExactMolWt_version
ExactMolWt.__doc__="""The exact molecular weight of the molecule

  >>> ExactMolWt(Chem.MolFromSmiles('CC'))
  30.04...
  >>> ExactMolWt(Chem.MolFromSmiles('[13CH3]C'))
  31.05...

"""


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

def NumRadicalElectrons(mol):
  """ The number of radical electrons the molecule has
    (says nothing about spin state)

  >>> NumRadicalElectrons(Chem.MolFromSmiles('CC'))
  0.0
  >>> NumRadicalElectrons(Chem.MolFromSmiles('C[CH3]'))
  0.0
  >>> NumRadicalElectrons(Chem.MolFromSmiles('C[CH2]'))
  1.0
  >>> NumRadicalElectrons(Chem.MolFromSmiles('C[CH]'))
  2.0
  >>> NumRadicalElectrons(Chem.MolFromSmiles('C[C]'))
  3.0
  

  """
  accum = 0.0
  for atom in mol.GetAtoms():
    accum += atom.GetNumRadicalElectrons()
  return accum
NumRadicalElectrons.version="1.0.0"

def _ChargeDescriptors(mol,force=False):
  if not force and hasattr(mol,'_chargeDescriptors'):
    return mol._chargeDescriptors
  chgs = rdPartialCharges.ComputeGasteigerCharges(mol)
  minChg=500.
  maxChg=-500.
  for at in mol.GetAtoms():
    chg = float(at.GetProp('_GasteigerCharge'))
    minChg = min(chg,minChg)
    maxChg = max(chg,maxChg)
  res = (minChg,maxChg)
  mol._chargeDescriptors=res
  return res
  

def MaxPartialCharge(mol,force=False):
  _,res = _ChargeDescriptors(mol,force)
  return res
MaxPartialCharge.version="1.0.0"

def MinPartialCharge(mol,force=False):
  res,_ = _ChargeDescriptors(mol,force)
  return res
MinPartialCharge.version="1.0.0"

def MaxAbsPartialCharge(mol,force=False):
  v1,v2 = _ChargeDescriptors(mol,force)
  return max(abs(v1),abs(v2))
MaxAbsPartialCharge.version="1.0.0"

def MinAbsPartialCharge(mol,force=False):
  v1,v2 = _ChargeDescriptors(mol,force)
  return min(abs(v1),abs(v2))
MinAbsPartialCharge.version="1.0.0"

from rdkit.Chem.EState.EState import MaxEStateIndex,MinEStateIndex,MaxAbsEStateIndex,MinAbsEStateIndex

_setupDescriptors(locals())



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
