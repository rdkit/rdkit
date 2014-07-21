## Automatically adapted for numpy.oldnumeric Jun 27, 2008 by -c

# $Id$
#
#  Copyright (C) 2000-2008  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" A module for molecules and stuff

 see Chem/index.html in the doc tree for documentation
 
"""
from rdkit import rdBase
from rdkit import RDConfig

from rdkit import DataStructs
from rdkit.Geometry import rdGeometry
from rdkit.Chem import PeriodicTable as pyPeriodicTable
from rdkit.Chem import rdchem
_HasSubstructMatchStr=rdchem._HasSubstructMatchStr
from rdkit.Chem.rdchem import *
from rdkit.Chem.rdmolfiles import *
from rdkit.Chem.rdmolops import *
from rdkit.Chem.inchi import *

def QuickSmartsMatch(smi,sma,unique=True,display=False):
  m = MolFromSmiles(smi)
  p = MolFromSmarts(sma)
  res = m.GetSubstructMatches(p,unique)
  if display:
    pass
  return res  

def CanonSmiles(smi,useChiral=1):
  m = MolFromSmiles(smi)
  return MolToSmiles(m,useChiral)

def SupplierFromFilename(fileN,delim='',**kwargs):
  ext = fileN.split('.')[-1].lower()
  if ext=='sdf':
    suppl = SDMolSupplier(fileN,**kwargs)
  elif ext=='csv':
    if not delim:
      delim = ','
    suppl = SmilesMolSupplier(fileN,delimiter=delim,**kwargs)
  elif ext=='txt':
    if not delim:
      delim='\t'
    suppl = SmilesMolSupplier(fileN,delimiter=delim,**kwargs)
  elif ext=='tdt':
    suppl = TDTMolSupplier(fileN,delimiter=delim,**kwargs)
  else:
    raise ValueError("unrecognized extension: %s"%ext)
    
  return suppl

def FindMolChiralCenters(mol,force=True,includeUnassigned=False):
  """
    >>> from rdkit import Chem
    >>> mol = Chem.MolFromSmiles('[C@H](Cl)(F)Br')
    >>> FindMolChiralCenters(mol)
    [(0, 'R')]
    >>> mol = Chem.MolFromSmiles('[C@@H](Cl)(F)Br')
    >>> FindMolChiralCenters(mol)
    [(0, 'S')]
  
    >>> FindMolChiralCenters(Chem.MolFromSmiles('CCC'))
    []

    By default unassigned stereo centers are not reported:
    >>> mol = Chem.MolFromSmiles('C[C@H](F)C(F)(Cl)Br')
    >>> FindMolChiralCenters(mol,force=True)
    [(1, 'S')]

    but this can be changed:
    >>> FindMolChiralCenters(mol,force=True,includeUnassigned=True)
    [(1, 'S'), (3, '?')]

    The handling of unassigned stereocenters for dependent stereochemistry is not correct:
    >>> Chem.FindMolChiralCenters(Chem.MolFromSmiles('C1CC(C)C(C)C(C)C1'),includeUnassigned=True)
    [(2, '?'), (6, '?')]
    >>> Chem.FindMolChiralCenters(Chem.MolFromSmiles('C1C[C@H](C)C(C)[C@H](C)C1'),includeUnassigned=True)
    [(2, 'S'), (4, '?'), (6, 'R')]
    
  """
  AssignStereochemistry(mol,force=force, flagPossibleStereoCenters=includeUnassigned)
  centers = []
  for atom in mol.GetAtoms():
    if atom.HasProp('_CIPCode'):
      centers.append((atom.GetIdx(),atom.GetProp('_CIPCode')))
    elif includeUnassigned and atom.HasProp('_ChiralityPossible'):
      centers.append((atom.GetIdx(),'?'))
  return centers

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
