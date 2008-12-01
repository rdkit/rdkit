## Automatically adapted for numpy.oldnumeric Jun 27, 2008 by -c

# $Id$
#
#  Copyright (C) 2006-2008  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
""" Import all RDKit chemistry modules
 
"""
import rdBase
import RDConfig
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
from rdSLNParse import *
import ForceField
Mol.Compute2DCoords = Compute2DCoords
Mol.ComputeGasteigerCharges = ComputeGasteigerCharges
import numpy

def TransformMol(mol,tform):
  """  Applies the transformation to a molecule and sets it up with
  a single conformer

  """
  newConf = Conformer()
  newConf.SetId(0)
  refConf = mol.GetConformer()
  for i in range(refConf.GetNumAtoms()):
    pos = list(refConf.GetAtomPosition(i))
    pos.append(1.0)
    newPos = numpy.dot(tform,numpy.array(pos))
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
  nAts = mol.GetNumAtoms()
  dm = []
  conf = mol.GetConformer(confId)
  for i in range(nAts):
    pi = conf.GetAtomPosition(i)
    for j in range(i+1,nAts):
      pj = conf.GetAtomPosition(j)
      dm.append((pi-pj).Length())
  dm = numpy.array(dm)
  apply(Compute2DCoordsMimicDistmat,(mol,dm),kwargs)
      
def GetBestRMS(ref,probe,refConfId=-1,probeConfId=-1,maps=None):
  if not maps:
    query = RemoveHs(probe)
    matches = ref.GetSubstructMatches(query)
    if not matches:
      raise ValueError,'mol %s does not match mol %s'%(ref.GetProp('_Name'),
                                                       probe.GetProp('_Name'))
    maps = []
    for match in matches:
      t=[]
      for j,idx in enumerate(match):
        t.append((j,idx))
    maps.append(t)
  bestRMS=1000.
  for amap in maps:
    rms=AlignMol(probe,ref,probeConfId,refConfId,atomMap=amap)
    if rms<bestRMS:
      bestRMS=rms
  return bestRMS

def EnumerateLibraryFromReaction(reaction,sidechainSets) :
  """ Returns a generator for the virtual library defined by
   a reaction and a sequence of sidechain sets

  >>> import Chem
  >>> from Chem import AllChem
  >>> s1=[Chem.MolFromSmiles(x) for x in ('NC','NCC')]
  >>> s2=[Chem.MolFromSmiles(x) for x in ('OC=O','OC(=O)C')]
  >>> rxn = AllChem.ReactionFromSmarts('[O:2]=[C:1][OH].[N:3]>>[O:2]=[C:1][N:3]')
  >>> r = AllChem.EnumerateLibraryFromReaction(rxn,[s2,s1])
  >>> [Chem.MolToSmiles(x[0]) for x in list(r)]
  ['CNC=O', 'CCNC=O', 'CNC(=O)C', 'CCNC(=O)C']

  Note that this is all done in a lazy manner, so "infinitely" large libraries can
  be done without worrying about running out of memory. Your patience will run out first:

  Define a set of 10000 amines:
  >>> amines = (Chem.MolFromSmiles('N'+'C'*x) for x in range(10000))

  ... a set of 10000 acids
  >>> acids = (Chem.MolFromSmiles('OC(=O)'+'C'*x) for x in range(10000))

  ... now the virtual library (1e8 compounds in principle):
  >>> r = AllChem.EnumerateLibraryFromReaction(rxn,[acids,amines])

  ... look at the first 4 compounds:
  >>> [Chem.MolToSmiles(r.next()[0]) for x in range(4)]
  ['NC=O', 'CNC=O', 'CCNC=O', 'CCCNC=O']

  
  """
  if len(sidechainSets) != reaction.GetNumReactantTemplates():
    raise ValueError,'%d sidechains provided, %d required'%(len(sidechainSets),reaction.getNumReactantTemplates())

  def _combiEnumerator(items,depth=0):
    for item in items[depth]:
      if depth+1 < len(items):
        v = _combiEnumerator(items,depth+1)
        for entry in v:
          l=[item]
          l.extend(entry)
          yield l
      else:
        yield [item]
  for chains in _combiEnumerator(sidechainSets):
    prodSets = reaction.RunReactants(chains)
    for prods in prodSets:
      yield prods



def ConstrainedEmbed(mol,core,useTethers,randomseed=2342):
  match = mol.GetSubstructMatch(core)
  if not match:
    raise ValueError,"molecule doesn't match the core"
  coordMap={}
  coreConf = core.GetConformer()
  for i,idxI in enumerate(match):
    corePtI = coreConf.GetAtomPosition(i)
    coordMap[idxI]=corePtI

  ci = EmbedMolecule(mol,coordMap=coordMap,randomSeed=randomseed)
  if ci<0:
    logger.error('could not embed molecule %s, no coordinates generated.'%mol.GetProp('_Name'))
  
  algMap=[]
  for i,itm in enumerate(match):
    algMap.append((itm,i))

  if not useTethers:
    # clean up the conformation
    ff = UFFGetMoleculeForceField(mol,confId=0)
    for i,idxI in enumerate(match):
      for j in range(i+1,len(match)):
        idxJ = match[j]
        d = coordMap[idxI].Distance(coordMap[idxJ])
        ff.AddDistanceConstraint(idxI,idxJ,d,d,100.)
    ff.Initialize()
    n=4
    more=ff.Minimize()
    while more and n:
      more=ff.Minimize()
      n-=1
    # rotate the embedded conformation onto the core:
    rms =AlignMol(mol,core,atomMap=algMap)
  else:
    # rotate the embedded conformation onto the core:
    rms = AlignMol(mol,core,atomMap=algMap)
    ff =  UFFGetMoleculeForceField(mol,confId=0)
    conf = core.GetConformer()
    for i in range(core.GetNumAtoms()):
      p =conf.GetAtomPosition(i)
      pIdx=ff.AddExtraPoint(p.x,p.y,p.z,fixed=True)-1
      ff.AddDistanceConstraint(pIdx,match[i],0,0,100.)
    ff.Initialize()
    n=4
    more=ff.Minimize(energyTol=1e-4,forceTol=1e-3)
    while more and n:
      more=ff.Minimize(energyTol=1e-4,forceTol=1e-3)
      n-=1
    # realign
    rms = AlignMol(mol,core,atomMap=algMap)
  print rms
  return mol


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
