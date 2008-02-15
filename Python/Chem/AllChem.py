# $Id$
#
#  Copyright (C) 2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
""" Import all RDKit chemistry modules
 
"""
import rdBase
import RDConfig
import Numeric
import DataStructs
from Geometry import rdGeometry
from Chem import *
from rdPartialCharges import *
from rdDepictor import *
from rdForceFieldHelpers import *
from Chem.ChemicalFeatures import *
from rdDistGeom import *
from rdMolAlign import *
from rdMolTransforms import *
from rdShapeHelpers import *
from rdChemReactions import *
import ForceField
Mol.Compute2DCoords = Compute2DCoords
Mol.ComputeGasteigerCharges = ComputeGasteigerCharges


def TransformMol(mol,tform):
  """  Applies the transformation to a molecule and sets it up with
  a single conformer

  """
  import Numeric
  newConf = Conformer()
  newConf.SetId(0)
  refConf = mol.GetConformer()
  for i in range(refConf.GetNumAtoms()):
    pos = list(refConf.GetAtomPosition(i))
    pos.append(1.0)
    newPos = Numeric.matrixmultiply(tform,Numeric.array(pos))
    newConf.SetAtomPosition(i,list(newPos)[:3])
  mol.RemoveAllConformers()
  mol.AddConformer(newConf)

def ComputeMolShape(mol,confId=-1,boxDim=(20,20,20),spacing=0.5,**kwargs):
  res = rdGeometry.UniformGrid3D(boxDim[0],boxDim[1],boxDim[2],spacing=spacing)
  apply(EncodeShape,(mol,res,confId),kwargs)
  return res
  
def ComputeMolVolume(mol,confId=-1,gridSpacing=0.1,boxMargin=2.0):
  import copy
  mol = copy.deepcopy(mol)
  conf = mol.GetConformer(confId)
  CanonicalizeConformer(conf)
  box = ComputeConfBox(conf)
  sideLen = ( box[1].x-box[0].x + 2*boxMargin, \
              box[1].y-box[0].y + 2*boxMargin, \
              box[1].z-box[0].z + 2*boxMargin )
  shape = rdGeometry.UniformGrid3D(sideLen[0],sideLen[1],sideLen[2],
                                   spacing=gridSpacing)
  EncodeShape(mol,shape,confId,ignoreHs=False,vdwScale=1.0)
  voxelVol = gridSpacing**3
  vol = 0.0
  occVect = shape.GetOccupancyVect()
  for i in range(len(occVect)):
    if occVect[i]==3: 
      vol+= voxelVol
  return vol

def GenerateDepictionMatching3DStructure(mol,reference,confId=-1,
                                         **kwargs):
  stripRef = RemoveHs(reference)
  nAts = mol.GetNumAtoms()
  dm = []
  conf = stripRef.GetConformer(confId)
  for i in range(nAts):
    pi = conf.GetAtomPosition(i)
    for j in range(i+1,nAts):
      pj = conf.GetAtomPosition(j)
      dm.append((pi-pj).Length())
  dm = Numeric.array(dm)
  apply(Compute2DCoordsMimicDistmat,(mol,dm),kwargs)
      
  
