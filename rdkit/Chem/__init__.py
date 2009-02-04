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

def GetSmartsMatchCDXML(mol,patt,maps,which=0,showAllAtoms=0):
  try:
    from rdkit.Chem import CDXMLWriter
  except:
    CDXMLWriter = None
  if CDXMLWriter is None:
    return ''
  try:
    from cStringIO import StringIO
  except:
    from StringIO import StringIO
  if not showAllAtoms and (which < 0 or not len(maps) or len(maps)<which+1):
    return ''
  io = StringIO()
  matchBonds = []
  matchAtoms = []
  if not showAllAtoms:
    mapL = maps[which]
    for bondIdx in range(patt.GetNumBonds()):
      bnd = patt.GetBondWithIdx(bondIdx)
      begIdx = bnd.GetBeginAtomIdx()
      endIdx = bnd.GetEndAtomIdx()
      molBegIdx = mapL[begIdx]
      molEndIdx = mapL[endIdx]
      matchBonds.append(mol.GetBondBetweenAtoms(molBegIdx,molEndIdx).GetIdx())
  else:
    for mapL in maps:
      for idx in mapL:
        if idx not in matchAtoms: matchAtoms.append(idx)

  CDXMLWriter.MolToCDXML(mol,io,highlightBonds=matchBonds,
                         highlightAtoms=matchAtoms)
  cdxml = io.getvalue()
  cdxml = cdxml.replace('><','>\n<')
  return cdxml

def DisplaySmartsMatch(mol,patt,maps,which=0,showAllAtoms=0):
  try:
    from rdkit.utils import chemdraw
  except:
    chemdraw = None
  cdxml = GetSmartsMatchCDXML(mol,patt,maps,which=which,showAllAtoms=showAllAtoms)
  if cdxml and chemdraw:
    chemdraw.CDXDisplay(cdxml)
  return cdxml
  
def QuickSmartsMatch(smi,sma,unique=1,display=0):
  m = MolFromSmiles(smi)
  p = MolFromSmarts(sma)
  res = m.GetSubstructMatches(p,unique)
  if display:
    DisplaySmartsMatch(m,p,res)
  return res  

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
  elif ext=='tdt':
    suppl = TDTMolSupplier(fileN,delimiter=delim,**kwargs)
  else:
    raise ValueError,"unrecognized extension: %s"%ext
    
  return suppl

def FindMolChiralCenters(mol):
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
  AssignAtomChiralCodes(mol)
  centers = []
  for atom in mol.GetAtoms():
    if atom.HasProp('_CIPCode'):
      centers.append((atom.GetIdx(),atom.GetProp('_CIPCode')))
  return centers
