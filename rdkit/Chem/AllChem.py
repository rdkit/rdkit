# $Id$
#
#  Copyright (C) 2006-2011  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" Import all RDKit chemistry modules
 
"""
from rdkit import rdBase
from rdkit import RDConfig
from rdkit import DataStructs
from rdkit.Geometry import rdGeometry
from rdkit.Chem import *
from rdkit.Chem.rdPartialCharges import *
from rdkit.Chem.rdDepictor import *
from rdkit.Chem.rdForceFieldHelpers import *
from rdkit.Chem.ChemicalFeatures import *
from rdkit.Chem.rdDistGeom import *
from rdkit.Chem.rdMolAlign import *
from rdkit.Chem.rdMolTransforms import *
from rdkit.Chem.rdShapeHelpers import *
from rdkit.Chem.rdChemReactions import *
try:
  from rdkit.Chem.rdSLNParse import *
except:
  pass
from rdkit.Chem.rdMolDescriptors import *
from rdkit import ForceField
Mol.Compute2DCoords = Compute2DCoords
Mol.ComputeGasteigerCharges = ComputeGasteigerCharges
import numpy

def TransformMol(mol,tform,confId=-1,keepConfs=False):
  """  Applies the transformation (usually a 4x4 double matrix) to a molecule
  if keepConfs is False then all but that conformer are removed

  """
  refConf = mol.GetConformer(confId)
  TransformConformer(refConf,tform)
  if not keepConfs:
    if confId==-1: confId=0
    allConfIds = [c.GetId() for c in mol.GetConformers()]
    for id in allConfIds:
      if not id==confId: mol.RemoveConformer(id)
    #reset the conf Id to zero since there is only one conformer left
    mol.GetConformer(confId).SetId(0) 

def ComputeMolShape(mol,confId=-1,boxDim=(20,20,20),spacing=0.5,**kwargs):
  """ returns a grid representation of the molecule's shape
  """
  res = rdGeometry.UniformGrid3D(boxDim[0],boxDim[1],boxDim[2],spacing=spacing)
  EncodeShape(mol,res,confId,**kwargs)
  return res
  
def ComputeMolVolume(mol,confId=-1,gridSpacing=0.2,boxMargin=2.0):
  """ Calculates the volume of a particular conformer of a molecule
  based on a grid-encoding of the molecular shape.

  """
  mol = rdchem.Mol(mol.ToBinary())
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
  occVect = shape.GetOccupancyVect()
  voxels = [1 for x in occVect if x==3]
  vol = voxelVol*len(voxels)
  return vol

def GenerateDepictionMatching2DStructure(mol,reference,confId=-1,
                                         referencePattern=None,
                                         acceptFailure=False,
                                         **kwargs):
  """ Generates a depiction for a molecule where a piece of the molecule
     is constrained to have the same coordinates as a reference.

     This is useful for, for example, generating depictions of SAR data
     sets so that the cores of the molecules are all oriented the same
     way.

  Arguments:
    - mol:          the molecule to be aligned, this will come back
                    with a single conformer.
    - reference:    a molecule with the reference atoms to align to;
                    this should have a depiction.
    - confId:       (optional) the id of the reference conformation to use
    - referencePattern:  (optional) an optional molecule to be used to 
                         generate the atom mapping between the molecule
                         and the reference.
    - acceptFailure: (optional) if True, standard depictions will be generated
                     for molecules that don't have a substructure match to the
                     reference; if False, a ValueError will be raised

  """
  if reference and referencePattern:
    if not reference.GetNumAtoms(onlyExplicit=True)==referencePattern.GetNumAtoms(onlyExplicit=True):
      raise ValueError,'When a pattern is provided, it must have the same number of atoms as the reference'
    referenceMatch = reference.GetSubstructMatch(referencePattern)
    if not referenceMatch:
      raise ValueError,"Reference does not map to itself"
  else:
    referenceMatch = range(reference.GetNumAtoms(onlyExplicit=True))
  if referencePattern:
    match = mol.GetSubstructMatch(referencePattern)
  else:
    match = mol.GetSubstructMatch(reference)

  if not match:
    if not acceptFailure:
      raise ValueError,'Substructure match with reference not found.'
    else:
      coordMap={}
  else:
    conf = reference.GetConformer()
    coordMap={}
    for i,idx in enumerate(match):
      pt3 = conf.GetAtomPosition(referenceMatch[i])
      pt2 = rdGeometry.Point2D(pt3.x,pt3.y)
      coordMap[idx] = pt2
  Compute2DCoords(mol,clearConfs=True,coordMap=coordMap,canonOrient=False)

def GenerateDepictionMatching3DStructure(mol,reference,confId=-1,
                                         **kwargs):
  """ Generates a depiction for a molecule where a piece of the molecule
     is constrained to have coordinates similar to those of a 3D reference
     structure.

  Arguments:
    - mol:          the molecule to be aligned, this will come back
                    with a single conformer.
    - reference:    a molecule with the reference atoms to align to;
                    this should have a depiction.
    - confId:       (optional) the id of the reference conformation to use

  """
  nAts = mol.GetNumAtoms()
  dm = []
  conf = reference.GetConformer(confId)
  for i in range(nAts):
    pi = conf.GetAtomPosition(i)
    #npi.z=0
    for j in range(i+1,nAts):
      pj = conf.GetAtomPosition(j)
      #pj.z=0
      dm.append((pi-pj).Length())
  dm = numpy.array(dm)
  Compute2DCoordsMimicDistmat(mol,dm,**kwargs)
      
def GetBestRMS(ref,probe,refConfId=-1,probeConfId=-1,maps=None):
  """ Returns the optimal RMS for aligning two molecules, taking
  symmetry into account. As a side-effect, the probe molecule is
  left in the aligned state.

  Arguments:
    - ref: the reference molecule
    - probe: the molecule to be aligned to the reference
    - refConfId: (optional) reference conformation to use
    - probeConfId: (optional) probe conformation to use
    - maps: (optional) a list of lists of (probeAtomId,refAtomId)
      tuples with the atom-atom mappings of the two molecules.
      If not provided, these will be generated using a substructure
      search.
  
  """
  if not maps:
    query = RemoveHs(probe)
    matches = ref.GetSubstructMatches(query,uniquify=False)
    if not matches:
      raise ValueError,'mol %s does not match mol %s'%(ref.GetProp('_Name'),
                                                       probe.GetProp('_Name'))
    maps = [list(enumerate(match)) for match in matches]
  bestRMS=1000.
  for amap in maps:
    rms=AlignMol(probe,ref,probeConfId,refConfId,atomMap=amap)
    if rms<bestRMS:
      bestRMS=rms
      bestMap = amap

  # finally repeate the best alignment :
  if bestMap != amap:
    AlignMol(probe,ref,probeConfId,refConfId,atomMap=bestMap)
  return bestRMS

def EnumerateLibraryFromReaction(reaction,sidechainSets) :
  """ Returns a generator for the virtual library defined by
   a reaction and a sequence of sidechain sets

  >>> from rdkit import Chem
  >>> from rdkit.Chem import AllChem
  >>> s1=[Chem.MolFromSmiles(x) for x in ('NC','NCC')]
  >>> s2=[Chem.MolFromSmiles(x) for x in ('OC=O','OC(=O)C')]
  >>> rxn = AllChem.ReactionFromSmarts('[O:2]=[C:1][OH].[N:3]>>[O:2]=[C:1][N:3]')
  >>> r = AllChem.EnumerateLibraryFromReaction(rxn,[s2,s1])
  >>> [Chem.MolToSmiles(x[0]) for x in list(r)]
  ['CNC=O', 'CCNC=O', 'CNC(C)=O', 'CCNC(C)=O']

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

def ConstrainedEmbed(mol,core,useTethers=True,coreConfId=-1,
                     randomseed=2342):
  """ generates an embedding of a molecule where part of the molecule
  is constrained to have particular coordinates

  Arguments
    - mol: the molecule to embed
    - core: the molecule to use as a source of constraints
    - useTethers: (optional) if True, the final conformation will be
        optimized subject to a series of extra forces that pull the
        matching atoms to the positions of the core atoms. Otherwise
        simple distance constraints based on the core atoms will be
        used in the optimization.
    - coreConfId: (optional) id of the core conformation to use
    - randomSeed: (optional) seed for the random number generator


    An example, start by generating a template with a 3D structure:
    >>> from rdkit.Chem import AllChem
    >>> template = AllChem.MolFromSmiles("c1nn(Cc2ccccc2)cc1")
    >>> AllChem.EmbedMolecule(template)
    0
    >>> AllChem.UFFOptimizeMolecule(template)
    0

    Here's a molecule:
    >>> mol = AllChem.MolFromSmiles("c1nn(Cc2ccccc2)cc1-c3ccccc3")

    Now do the constrained embedding
    >>> newmol=AllChem.ConstrainedEmbed(mol, template)

    Demonstrate that the positions are the same:
    >>> newp=newmol.GetConformer().GetAtomPosition(0)
    >>> molp=mol.GetConformer().GetAtomPosition(0)
    >>> list(newp-molp)==[0.0,0.0,0.0]
    True
    >>> newp=newmol.GetConformer().GetAtomPosition(1)
    >>> molp=mol.GetConformer().GetAtomPosition(1)
    >>> list(newp-molp)==[0.0,0.0,0.0]
    True

  """
  match = mol.GetSubstructMatch(core)
  if not match:
    raise ValueError,"molecule doesn't match the core"
  coordMap={}
  coreConf = core.GetConformer(coreConfId)
  for i,idxI in enumerate(match):
    corePtI = coreConf.GetAtomPosition(i)
    coordMap[idxI]=corePtI

  ci = EmbedMolecule(mol,coordMap=coordMap,randomSeed=randomseed)
  if ci<0:
    raise ValueError,'Could not embed molecule.'
  
  algMap=[(j,i) for i,j in enumerate(match)]

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
  mol.SetProp('EmbedRMS',str(rms))
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
