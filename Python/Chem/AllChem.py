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
      
  
