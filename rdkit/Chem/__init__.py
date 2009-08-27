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

from rdkit.Geometry import rdGeometry
import rdchem
_HasSubstructMatchStr=rdchem._HasSubstructMatchStr
from rdchem import *
from rdmolfiles import *
from rdmolops import *

def CanonSmiles(smi,useChiral=1):
  m = MolFromSmiles(smi)
  return MolToSmiles(m,useChiral)

def SmilesRoundtrip(smi,useChiral=1):
  m = MolFromSmiles(smi)
  refSmi = MolToSmiles(m,useChiral)
  m2 = MolFromSmiles(refSmi)
  smi = MolToSmiles(m2,useChiral)
  if smi!=refSmi:
    print refSmi
    print smi
  return refSmi==smi

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
