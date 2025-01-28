#
#  Copyright (C) 2006-2022  greg Landrum and other RDKit contributors
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" Import all RDKit chemistry modules

"""
import sys
import warnings
from collections import namedtuple

import numpy

from rdkit import DataStructs, ForceField, RDConfig, rdBase
from rdkit.Chem import *
from rdkit.Chem.ChemicalFeatures import *
from rdkit.Chem.EnumerateStereoisomers import (EnumerateStereoisomers, StereoEnumerationOptions)
from rdkit.Chem.rdChemReactions import *
from rdkit.Chem.rdDepictor import *
from rdkit.Chem.rdDistGeom import *
from rdkit.Chem.rdFingerprintGenerator import *
from rdkit.Chem.rdForceFieldHelpers import *
from rdkit.Chem.rdMolAlign import *
from rdkit.Chem.rdMolDescriptors import *
from rdkit.Chem.rdMolEnumerator import *
from rdkit.Chem.rdMolTransforms import *
from rdkit.Chem.rdPartialCharges import *
from rdkit.Chem.rdqueries import *
from rdkit.Chem.rdReducedGraphs import *
from rdkit.Chem.rdShapeHelpers import *
from rdkit.Geometry import rdGeometry
from rdkit.RDLogger import logger

try:
  from rdkit.Chem.rdSLNParse import *
except ImportError:
  pass

Mol.Compute2DCoords = Compute2DCoords
Mol.ComputeGasteigerCharges = ComputeGasteigerCharges
logger = logger()


def TransformMol(mol, tform, confId=-1, keepConfs=False):
  """  Applies the transformation (usually a 4x4 double matrix) to a molecule
    
  Arguments:
    - mol: the molecule to be transformed
    - tform: the transformation to apply
    - confId: (optional) the conformer id to transform
    - keepConfs: (optional) if keepConfs is False then all but that conformer are removed
  
  """
  refConf = mol.GetConformer(confId)
  TransformConformer(refConf, tform)
  if not keepConfs:
    if confId == -1:
      confId = 0
    allConfIds = [c.GetId() for c in mol.GetConformers()]
    for cid in allConfIds:
      if not cid == confId:
        mol.RemoveConformer(cid)
    # reset the conf Id to zero since there is only one conformer left
    mol.GetConformer(confId).SetId(0)


def ComputeMolShape(mol, confId=-1, boxDim=(20, 20, 20), spacing=0.5, **kwargs):
  """ returns a grid representation of the molecule's shape

  Arguments:
    - mol: the molecule
    - confId: (optional) the conformer id to use
    - boxDim: (optional) the dimensions of the box to use
    - spacing: (optional) the spacing to use
    - kwargs: additional arguments to pass to the encoding function

  Returns:
    a UniformGrid3D object
  """
  res = rdGeometry.UniformGrid3D(boxDim[0], boxDim[1], boxDim[2], spacing=spacing)
  EncodeShape(mol, res, confId, **kwargs)
  return res


def ComputeMolVolume(mol, confId=-1, gridSpacing=0.2, boxMargin=2.0):
  """ Calculates the volume of a particular conformer of a molecule
    based on a grid-encoding of the molecular shape.


  Arguments:
    - mol: the molecule
    - confId: (optional) the conformer id to use
    - gridSpacing: (optional) the spacing to use 
    - boxMargin: (optional) the margin to use around the molecule

    A bit of demo as well as a test of github #1883:

    >>> from rdkit import Chem
    >>> from rdkit.Chem import AllChem
    >>> mol = Chem.AddHs(Chem.MolFromSmiles('C'))
    >>> AllChem.EmbedMolecule(mol)
    0
    >>> ComputeMolVolume(mol)
    28...
    >>> mol = Chem.AddHs(Chem.MolFromSmiles('O'))
    >>> AllChem.EmbedMolecule(mol)
    0
    >>> ComputeMolVolume(mol)
    20...

    """
  mol = rdchem.Mol(mol)
  conf = mol.GetConformer(confId)
  CanonicalizeConformer(conf, ignoreHs=False)
  box = ComputeConfBox(conf)
  sideLen = (box[1].x - box[0].x + 2 * boxMargin, box[1].y - box[0].y + 2 * boxMargin,
             box[1].z - box[0].z + 2 * boxMargin)
  shape = rdGeometry.UniformGrid3D(sideLen[0], sideLen[1], sideLen[2], spacing=gridSpacing)
  EncodeShape(mol, shape, confId, ignoreHs=False, vdwScale=1.0)
  voxelVol = gridSpacing**3
  occVect = shape.GetOccupancyVect()
  voxels = [1 for x in occVect if x == 3]
  vol = voxelVol * len(voxels)
  return vol


def GetConformerRMS(mol, confId1, confId2, atomIds=None, prealigned=False):
  """ Returns the RMS between two conformations.
    By default, the conformers will be aligned to the first conformer
    before the RMS calculation and, as a side-effect, the second will be left
    in the aligned state.

    Arguments:
      - mol:        the molecule
      - confId1:    the id of the first conformer
      - confId2:    the id of the second conformer
      - atomIds:    (optional) list of atom ids to use a points for
                    alingment - defaults to all atoms
      - prealigned: (optional) by default the conformers are assumed
                    be unaligned and the second conformer be aligned
                    to the first

    """
  # align the conformers if necessary
  # Note: the reference conformer is always the first one
  if not prealigned:
    if atomIds:
      AlignMolConformers(mol, confIds=[confId1, confId2], atomIds=atomIds)
    else:
      AlignMolConformers(mol, confIds=[confId1, confId2])

  # calculate the RMS between the two conformations
  conf1 = mol.GetConformer(id=confId1)
  conf2 = mol.GetConformer(id=confId2)
  ssr = 0
  for i in range(mol.GetNumAtoms()):
    d = conf1.GetAtomPosition(i).Distance(conf2.GetAtomPosition(i))
    ssr += d * d
  ssr /= mol.GetNumAtoms()
  return numpy.sqrt(ssr)


def GetConformerRMSMatrix(mol, atomIds=None, prealigned=False):
  """ Returns the RMS matrix of the conformers of a molecule.
    As a side-effect, the conformers will be aligned to the first
    conformer (i.e. the reference) and will left in the aligned state.

    Arguments:
      - mol:     the molecule
      - atomIds: (optional) list of atom ids to use a points for
                 alingment - defaults to all atoms
      - prealigned: (optional) by default the conformers are assumed
                    be unaligned and will therefore be aligned to the
                    first conformer

    Note that the returned RMS matrix is symmetrical, i.e. it is the
    lower half of the matrix, e.g. for 5 conformers::

      rmsmatrix = [ a,
                    b, c,
                    d, e, f,
                    g, h, i, j]

    where a is the RMS between conformers 0 and 1, b is the RMS between
    conformers 0 and 2, etc.
    This way it can be directly used as distance matrix in e.g. Butina
    clustering.

    """
  # if necessary, align the conformers
  # Note: the reference conformer is always the first one
  rmsvals = []
  confIds = [conf.GetId() for conf in mol.GetConformers()]
  if not prealigned:
    if atomIds:
      AlignMolConformers(mol, atomIds=atomIds, RMSlist=rmsvals)
    else:
      AlignMolConformers(mol, RMSlist=rmsvals)
  else:  # already prealigned
    for i in range(1, len(confIds)):
      rmsvals.append(
        GetConformerRMS(mol, confIds[0], confIds[i], atomIds=atomIds, prealigned=prealigned))
  # loop over the conformations (except the reference one)
  cmat = []
  for i in range(1, len(confIds)):
    cmat.append(rmsvals[i - 1])
    for j in range(1, i):
      cmat.append(GetConformerRMS(mol, confIds[i], confIds[j], atomIds=atomIds, prealigned=True))
  return cmat


def EnumerateLibraryFromReaction(reaction, sidechainSets, returnReactants=False):
  """ Returns a generator for the virtual library defined by
    a reaction and a sequence of sidechain sets

    Arguments:
      - reaction: the reaction
      - sidechainSets: a sequence of sequences of sidechains
      - returnReactants: (optional) if True, the generator will return information about the reactants
                         as well as the products

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

    >>> [Chem.MolToSmiles(next(r)[0]) for x in range(4)]
    ['NC=O', 'CNC=O', 'CCNC=O', 'CCCNC=O']

    Here's what returnReactants does:

    >>> l = list(AllChem.EnumerateLibraryFromReaction(rxn,[s2,s1],returnReactants=True))
    >>> type(l[0])
    <class 'rdkit.Chem.AllChem.ProductReactants'>
    >>> [Chem.MolToSmiles(x) for x in l[0].reactants]
    ['O=CO', 'CN']
    >>> [Chem.MolToSmiles(x) for x in l[0].products]
    ['CNC=O']

    """
  if len(sidechainSets) != reaction.GetNumReactantTemplates():
    raise ValueError('%d sidechains provided, %d required' %
                     (len(sidechainSets), reaction.GetNumReactantTemplates()))

  def _combiEnumerator(items, depth=0):
    for item in items[depth]:
      if depth + 1 < len(items):
        v = _combiEnumerator(items, depth + 1)
        for entry in v:
          l = [item]
          l.extend(entry)
          yield l
      else:
        yield [item]

  ProductReactants = namedtuple('ProductReactants', 'products,reactants')
  for chains in _combiEnumerator(sidechainSets):
    prodSets = reaction.RunReactants(chains)
    for prods in prodSets:
      if returnReactants:
        yield ProductReactants(prods, chains)
      else:
        yield prods


def ConstrainedEmbed(mol, core, useTethers=True, coreConfId=-1, randomseed=2342,
                     getForceField=UFFGetMoleculeForceField, **kwargs):
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
      - getForceField: (optional) a function to use to get a force field
          for the final cleanup
      - kwargs: additional arguments to pass to the embedding function

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
  
    >>> mol = AllChem.ConstrainedEmbed(mol, template)

    Demonstrate that the positions are nearly the same with template:

    >>> import math
    >>> molp = mol.GetConformer().GetAtomPosition(0)
    >>> templatep = template.GetConformer().GetAtomPosition(0)
    >>> all(math.isclose(v, 0.0, abs_tol=0.01) for v in molp-templatep)
    True
    >>> molp = mol.GetConformer().GetAtomPosition(1)
    >>> templatep = template.GetConformer().GetAtomPosition(1)
    >>> all(math.isclose(v, 0.0, abs_tol=0.01) for v in molp-templatep)
    True

    """
  match = mol.GetSubstructMatch(core)
  if not match:
    raise ValueError("molecule doesn't match the core")
  coordMap = {}
  coreConf = core.GetConformer(coreConfId)
  for i, idxI in enumerate(match):
    corePtI = coreConf.GetAtomPosition(i)
    coordMap[idxI] = corePtI

  ci = EmbedMolecule(mol, coordMap=coordMap, randomSeed=randomseed, **kwargs)
  if ci < 0:
    raise ValueError('Could not embed molecule.')

  algMap = [(j, i) for i, j in enumerate(match)]

  if not useTethers:
    # clean up the conformation
    ff = getForceField(mol, confId=0)
    for i, idxI in enumerate(match):
      for j in range(i + 1, len(match)):
        idxJ = match[j]
        d = coordMap[idxI].Distance(coordMap[idxJ])
        ff.AddDistanceConstraint(idxI, idxJ, d, d, 100.)
    ff.Initialize()
    n = 4
    more = ff.Minimize()
    while more and n:
      more = ff.Minimize()
      n -= 1
    # rotate the embedded conformation onto the core:
    rms = AlignMol(mol, core, atomMap=algMap)
  else:
    # rotate the embedded conformation onto the core:
    rms = AlignMol(mol, core, atomMap=algMap)
    ff = getForceField(mol, confId=0)
    conf = core.GetConformer()
    for i in range(core.GetNumAtoms()):
      p = conf.GetAtomPosition(i)
      pIdx = ff.AddExtraPoint(p.x, p.y, p.z, fixed=True) - 1
      ff.AddDistanceConstraint(pIdx, match[i], 0, 0, 100.)
    ff.Initialize()
    n = 4
    more = ff.Minimize(energyTol=1e-4, forceTol=1e-3)
    while more and n:
      more = ff.Minimize(energyTol=1e-4, forceTol=1e-3)
      n -= 1
    # realign
    rms = AlignMol(mol, core, atomMap=algMap)
  mol.SetProp('EmbedRMS', str(rms))
  return mol


def AssignBondOrdersFromTemplate(refmol, mol):
  """ assigns bond orders to a molecule based on the
    bond orders in a template molecule

    Arguments
      - refmol: the template molecule
      - mol: the molecule to assign bond orders to

    An example, start by generating a template from a SMILES
    and read in the PDB structure of the molecule

    >>> import os
    >>> from rdkit.Chem import AllChem
    >>> template = AllChem.MolFromSmiles("CN1C(=NC(C1=O)(c2ccccc2)c3ccccc3)N")
    >>> mol = AllChem.MolFromPDBFile(os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data', '4DJU_lig.pdb'))
    >>> len([1 for b in template.GetBonds() if b.GetBondTypeAsDouble() == 1.0])
    8
    >>> len([1 for b in mol.GetBonds() if b.GetBondTypeAsDouble() == 1.0])
    22

    Now assign the bond orders based on the template molecule

    >>> newMol = AllChem.AssignBondOrdersFromTemplate(template, mol)
    >>> len([1 for b in newMol.GetBonds() if b.GetBondTypeAsDouble() == 1.0])
    8

    Note that the template molecule should have no explicit hydrogens
    else the algorithm will fail.

    It also works if there are different formal charges (this was github issue 235):

    >>> template=AllChem.MolFromSmiles('CN(C)C(=O)Cc1ccc2c(c1)NC(=O)c3ccc(cc3N2)c4ccc(c(c4)OC)[N+](=O)[O-]')
    >>> mol = AllChem.MolFromMolFile(os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data', '4FTR_lig.mol'))
    >>> AllChem.MolToSmiles(mol)
    'COC1CC(C2CCC3C(O)NC4CC(CC(O)N(C)C)CCC4NC3C2)CCC1N(O)O'
    >>> newMol = AllChem.AssignBondOrdersFromTemplate(template, mol)
    >>> AllChem.MolToSmiles(newMol)
    'COc1cc(-c2ccc3c(c2)Nc2ccc(CC(=O)N(C)C)cc2NC3=O)ccc1[N+](=O)[O-]'

    """
  refmol2 = rdchem.Mol(refmol)
  mol2 = rdchem.Mol(mol)
  # do the molecules match already?
  matching = mol2.GetSubstructMatch(refmol2)
  if not matching:  # no, they don't match
    # check if bonds of mol are SINGLE
    for b in mol2.GetBonds():
      if b.GetBondType() != BondType.SINGLE:
        b.SetBondType(BondType.SINGLE)
        b.SetIsAromatic(False)
    # set the bonds of mol to SINGLE
    for b in refmol2.GetBonds():
      b.SetBondType(BondType.SINGLE)
      b.SetIsAromatic(False)
    # set atom charges to zero;
    for a in refmol2.GetAtoms():
      a.SetFormalCharge(0)
    for a in mol2.GetAtoms():
      a.SetFormalCharge(0)

    matching = mol2.GetSubstructMatches(refmol2, uniquify=False)
    # do the molecules match now?
    if matching:
      if len(matching) > 1:
        logger.warning("More than one matching pattern found - picking one")
      matching = matching[0]
      # apply matching: set bond properties
      for b in refmol.GetBonds():
        atom1 = matching[b.GetBeginAtomIdx()]
        atom2 = matching[b.GetEndAtomIdx()]
        b2 = mol2.GetBondBetweenAtoms(atom1, atom2)
        b2.SetBondType(b.GetBondType())
        b2.SetIsAromatic(b.GetIsAromatic())
      # apply matching: set atom properties
      for a in refmol.GetAtoms():
        a2 = mol2.GetAtomWithIdx(matching[a.GetIdx()])
        a2.SetHybridization(a.GetHybridization())
        a2.SetIsAromatic(a.GetIsAromatic())
        a2.SetNumExplicitHs(a.GetNumExplicitHs())
        a2.SetFormalCharge(a.GetFormalCharge())
      SanitizeMol(mol2)
      if hasattr(mol2, '__sssAtoms'):
        mol2.__sssAtoms = None  # we don't want all bonds highlighted
    else:
      raise ValueError("No matching found")
  return mol2
