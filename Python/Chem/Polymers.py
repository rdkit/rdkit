# $Id$
#
# Copyright (C) 2002-2006 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
""" Python functions for manipulating polymers

"""
raise NotImplementedError
from Chem import Graphs,Lipinski
import exceptions

class PolymerError(ValueError):
  pass

_linkerDataName='polymer_linkers'
_backboneDataName='polymer_backbone'

#---------------------------------------------------------------------------------
#
#   Infrastructure
#
#---------------------------------------------------------------------------------

def FindLinkers(mol,linkMarkers='X',linkerDataName=_linkerDataName):
  """ *Internal use only*
    *NOTES:*

      - the molecule is modified

  """
  if mol.HasData(linkerDataName):
    linkers = eval(mol.GetData(linkerDataName))
  else:
    linkers = []
    for atom in mol.GetAtoms():
      if atom.GetSymbol() == linkMarkers:
        linkers.append(atom.GetIdx())
    mol.SetData(linkerDataName,repr(linkers))
  return linkers

def FindBackboneAtoms(mol):
  """  *Internal use only*

    *NOTES:*

      - if the mol already has linker atoms, we'll use those,
        otherwise we're gonna calculate them

      - the molecule is modified

  """
  if mol.HasData(_backboneDataName):
    linkPath = eval(mol.GetData(_backboneDataName))
  else:    
    linkers = FindLinkers(mol)
    if len(linkers) != 2:
      raise PolymerError,'monomers must have two linkers'

    dMat,pathMat = mol.ShortestPaths()
    linkPath = Graphs.GetPath(pathMat,linkers[0]-1,linkers[1]-1)
    linkPath = map(lambda x:x+1,linkPath)

    mol.SetData(_backboneDataName,repr(linkPath))
  return linkPath


#---------------------------------------------------------------------------------
#
#   Descriptors
#
#---------------------------------------------------------------------------------
def FindAtomsInBackbone(mol,smarts,backbone=None,acceptPartial=0,removeLinkers=1):
  """ finds the atoms in a polymer backbone that match a SMARTS pattern

  *ARGUMENTS:*

    - mol: a Mol instance

    - smarts: a SmartsPattern instance

    - backbone: (optional) a sequence with the integer indices of atoms
       in the polymer backbone

    - acceptPartial: (optional) if this is true then multi-atom matches which include
       some (but not all) backbone atoms are accepted

    - removeLinkers: (optional) if this is not set, then matches which include the linker
       atoms are also possible

  *RETURNS:*

      a list of lists with matches

  """
  res = []
  if backbone is None:
    backbone = FindBackboneAtoms(mol)
    
  if removeLinkers:
    backbone = backbone[:]
    for linker in FindLinkers(mol):
      backbone.remove(linker)

  smarts.Match(mol,0)
  for match in smarts.GetUMapList():
    keep = len(match)
    for idx in match:
      if idx not in backbone:
        keep -= 1
        if not acceptPartial:
          break
    if keep == len(match) or (acceptPartial and keep):  
      res.append(match)
  return res    

def FindAtomsInSidechain(mol,smarts,backbone=None,acceptPartial=0):
  """ finds the atoms in a polymer sidechain that match a SMARTS pattern

  *ARGUMENTS:*

    - mol: a Mol instance

    - smarts: a SmartsPattern instance

    - backbone: (optional) a sequence with the integer indices of atoms
       in the polymer backbone

    - acceptPartial: (optional) if this is true then multi-atom matches which include
       some (but not all) backbone atoms are accepted

  *RETURNS:*

      a list of lists with matches

  """
  res = []
  if backbone is None:
    backbone = FindBackboneAtoms(mol)

  smarts.Match(mol,0)
  for match in smarts.GetUMapList():
    keep = len(match)
    for idx in match:
      if idx in backbone:
        keep -= 1
        if not acceptPartial:
          break
    if keep == len(match) or (acceptPartial and keep):  
      res.append(match)
  return res    

def FindRotorsInBackbone(mol,backbone=None):
  """ self-explanatory, uses _Chem.Lipinski_ definition
  """
  return FindAtomsInBackbone(mol,Lipinski.RotatableBondSmarts,backbone)
def FindRotorsInSidechain(mol,backbone=None):
  """ self-explanatory, uses _Chem.Lipinski_ definition
  """
  return FindAtomsInSidechain(mol,Lipinski.RotatableBondSmarts,backbone,acceptPartial=1)

def FindHeteroatomsInBackbone(mol,backbone=None):
  """ self-explanatory, uses _Chem.Lipinski_ definition
  """
  return FindAtomsInBackbone(mol,Lipinski.HeteroatomSmarts,backbone)
def FindHeteroatomsInSidechain(mol,backbone=None):
  """ self-explanatory, uses _Chem.Lipinski_ definition
  """
  return FindAtomsInSidechain(mol,Lipinski.HeteroatomSmarts,backbone)

def FindHDonorsInBackbone(mol,backbone=None):
  """ self-explanatory, uses _Chem.Lipinski_ definition
  """
  return FindAtomsInBackbone(mol,Lipinski.HDonorSmarts,backbone)
def FindHDonorsInSidechain(mol,backbone=None):
  """ self-explanatory, uses _Chem.Lipinski_ definition
  """
  return FindAtomsInSidechain(mol,Lipinski.HDonorSmarts,backbone)

def FindHAcceptorsInBackbone(mol,backbone=None):
  """ self-explanatory, uses _Chem.Lipinski_ definition
  """
  return FindAtomsInBackbone(mol,Lipinski.HAcceptorSmarts,backbone)
def FindHAcceptorsInSidechain(mol,backbone=None):
  """ self-explanatory, uses _Chem.Lipinski_ definition
  """
  return FindAtomsInSidechain(mol,Lipinski.HAcceptorSmarts,backbone)


