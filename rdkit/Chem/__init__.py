## Automatically adapted for numpy.oldnumeric Jun 27, 2008 by -c

# $Id$
#
#  Copyright (C) 2000-2008  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
""" A module for molecules and stuff

 see Chem/index.html in the doc tree for documentation
 
"""
from rdkit import rdBase
from rdkit import RDConfig

from rdkit import DataStructs
from rdkit.Geometry import rdGeometry
import PeriodicTable as pyPeriodicTable
import rdchem
_HasSubstructMatchStr=rdchem._HasSubstructMatchStr
from rdchem import *
from rdmolfiles import *
from rdmolops import *

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
    raise ValueError,"unrecognized extension: %s"%ext
    
  return suppl

def FindMolChiralCenters(mol,force=True):
  """
    >>> mol = Chem.MolFromSmiles('[C@H](Cl)(F)Br')
    >>> FindMolChiralCenters(mol)
    [(0, 'R')]
    >>> mol = Chem.MolFromSmiles('[C@@H](Cl)(F)Br')
    >>> FindMolChiralCenters(mol)
    [(0, 'S')]
  
    >>> FindMolChiralCenters(Chem.MolFromSmiles('CCC'))
    []
  
  """
  AssignStereochemistry(mol,force=force)
  centers = []
  for atom in mol.GetAtoms():
    if atom.HasProp('_CIPCode'):
      centers.append((atom.GetIdx(),atom.GetProp('_CIPCode')))
  return centers
