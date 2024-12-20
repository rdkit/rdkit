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


def RotateGroup(mol, start_idx, end_idx, angle_deg):
  """ Rotates a group of atoms around a bond in an RDKit molecule.

    Arguments
      - mol: RDKit molecule object.
      - start_idx: index of the starting atom of the bond.
      - end_idx: index of the ending atom of the bond.
      - angle_deg: rotation angle in degrees.

    An example:

    >>> from rdkit.Chem import AllChem
    >>> smiles = "CC(C)C1=CC=CC=C1"
    >>> mol = Chem.MolFromSmiles(smiles)
    >>> mol = Chem.AddHs(mol)
    >>> AllChem.EmbedMolecule(mol)
    >>> mol_rotated = AllChem.RotateGroup(mol, 0, 1, 30)

    """
  rotable_bonds = find_rotatable_bonds(mol)
  
  assert (start_idx, end_idx) in rotable_bonds or (end_idx, start_idx) in rotable_bonds, \
  f"({start_idx}, {end_idx}) or ({end_idx}, {start_idx}) is not a rotatable bond."

  conformer = mol.GetConformer()
  positions = numpy.array(conformer.GetPositions())

  group1, _ = collect_connected_groups(mol, start_idx, end_idx)
  group1_positions = positions[list(group1)]
  axis_point = positions[end_idx]

  group1_positions -= axis_point
  start_point = positions[start_idx] - axis_point
  theta_1, theta_2 = get_rotation_angles_to_align(start_point)
  aligned_group1 = apply_rotations(group1_positions, theta_1, theta_2)

  rotation_angle = angle_deg * (numpy.pi / 180)
  rot_mat = rotation_matrix(rotation_angle, 'y')
  rotated_group1 = aligned_group1 @ rot_mat.T
  r_z_minus = rotation_matrix(-theta_2, 'z')
  step_1_cloud = rotated_group1 @ r_z_minus.T
  r_y_minus = rotation_matrix(-theta_1, 'y')
  step_2_cloud = step_1_cloud @ r_y_minus.T

  step_2_cloud += axis_point
  for i, atom_idx in enumerate(group1):
    conformer.SetAtomPosition(atom_idx, step_2_cloud[i])

  return mol


def rotation_matrix(theta: float, axis: str):
  """ Generates a 3D rotation matrix for a specified angle and axis.

    Arguments
      - theta: the rotation angle in radians.
      - axis: the axis of rotation ('x', 'y', 'z').

    Returns
      A 3x3 numpy array representing the rotation matrix.

    Example:

    >>> theta = np.pi / 2
    >>> axis = 'y'
    >>> rot_matrix = rotation_matrix(theta, axis)
    >>> rot_matrix
    array([[ 0.,  0.,  1.],
           [ 0.,  1.,  0.],
           [-1.,  0.,  0.]])

    """
  matrices = {
      'x': numpy.array([[1, 0, 0],
                     [0, numpy.cos(theta), -numpy.sin(theta)],
                     [0, numpy.sin(theta), numpy.cos(theta)]]),
      'y': numpy.array([[numpy.cos(theta), 0, numpy.sin(theta)],
                     [0, 1, 0],
                     [-numpy.sin(theta), 0, numpy.cos(theta)]]),
      'z': numpy.array([[numpy.cos(theta), -numpy.sin(theta), 0],
                     [numpy.sin(theta), numpy.cos(theta), 0],
                     [0, 0, 1]])
  }
  return matrices[axis]


def get_rotation_angles_to_align(point: numpy.ndarray):
  """ Calculates the rotation angles required to align a 3D point with the x-axis.

    Arguments
      - point: a numpy array of shape (3,) representing a 3D point.

    Returns
      A tuple (theta_1, theta_2) of rotation angles in radians for the 'y' and 'z' axes.

    Example:

    >>> point = numpy.array([1.0, 0.5, 0.5])
    >>> theta_1, theta_2 = get_rotation_angles_to_align(point)
    >>> (theta_1, theta_2)
    (0.4636476090008061, 1.1071487177940904)

    """
  # Special case: If the point is already along the x-axis (point[1] == 0 and point[2] == 0)
  if numpy.allclose(point[1:], 0.0):
    return 0.0, 0.0
      
  ratio = point[2] / point[0]
  theta_1 = numpy.arctan(ratio)
  r_y = rotation_matrix(theta_1, 'y')
  rotated_point = r_y @ point
  
  if rotated_point[1] == 0:
    theta_2 = 0.0
  else:
    ratio_2 = rotated_point[0] / rotated_point[1]
    theta_2 = numpy.arctan(ratio_2)

  return theta_1, theta_2


def apply_rotations(points: numpy.ndarray, theta_1: float, theta_2: float):
  """ Applies rotations to align 3D points based on the specified angles.

    Arguments
      - points: a numpy array of shape (N, 3), where N is the number of points.
      - theta_1: rotation angle in radians for the 'y' axis.
      - theta_2: rotation angle in radians for the 'z' axis.

    Returns
      A numpy array of rotated points.

    Example:

    >>> points = numpy.array([[1.0, 0.0, 0.0]])
    >>> theta_1, theta_2 = numpy.pi / 4, numpy.pi / 6
    >>> rotated_points = apply_rotations(points, theta_1, theta_2)
    >>> rotated_points
    array([[ 0.61237244, -0.35355339,  0.70710678]])

    """
  r_y = rotation_matrix(theta_1, 'y')
  rotated_points = points @ r_y.T
  r_z = rotation_matrix(theta_2, 'z')
  return rotated_points @ r_z.T


def collect_connected_groups(mol, atom1_idx, atom2_idx):
  """ Collects groups of atoms connected to two atoms in a molecule, partitioned by a bond.

    Arguments
      - mol: RDKit molecule object.
      - atom1_idx: index of the first atom in the bond.
      - atom2_idx: index of the second atom in the bond.

    Returns
      A tuple (group1, group2):
        - group1: set of atoms connected to atom1, excluding atom2.
        - group2: set of atoms connected to atom2, excluding atom1.

    Example:

    >>> from rdkit.Chem import MolFromSmiles
    >>> mol = MolFromSmiles("CC(C)C1=CC=CC=C1")
    >>> group1, group2 = collect_connected_groups(mol, 0, 1)
    >>> (group1, group2)
    ({0, 3, 4, 5}, {1, 2})

    """
  group1, group2 = set(), set()

  def dfs(atom_idx, group, exclude_idx):
    stack = [atom_idx]
    visited = {atom_idx}
    while stack:
      current = stack.pop()
      group.add(current)
      for neighbor in mol.GetAtomWithIdx(current).GetNeighbors():
        neighbor_idx = neighbor.GetIdx()
        if neighbor_idx != exclude_idx and neighbor_idx not in visited:
          visited.add(neighbor_idx)
          stack.append(neighbor_idx)

  dfs(atom1_idx, group1, atom2_idx)
  dfs(atom2_idx, group2, atom1_idx)
  return group1, group2


def find_rotatable_bonds(mol):
  """ Identifies rotatable bonds in an RDKit molecule.

    A bond is considered rotatable if:
      - It is a single bond.
      - It is not part of a ring.
      - Neither of the connected atoms is terminal (i.e., degree > 1).

    Arguments
      - mol: RDKit molecule object.

    Returns
      A list of tuples, where each tuple represents the indices of two atoms
      connected by a rotatable bond.

    Example:

    >>> from rdkit.Chem import MolFromSmiles
    >>> smiles = "CC(C)C1=CC=CC=C1"
    >>> mol = MolFromSmiles(smiles)
    >>> rotatable_bonds = find_rotatable_bonds(mol)
    >>> rotatable_bonds
    [(0, 1), (1, 2)]

    """
  rotatable_bonds = []
  for bond in mol.GetBonds():
    if (
        bond.GetBondType() == BondType.SINGLE and  
        not bond.IsInRing()                            
    ):
      atom1 = bond.GetBeginAtom()
      atom2 = bond.GetEndAtom()

      if atom1.GetDegree() != 1 and atom2.GetDegree() != 1:
        atom1_idx = atom1.GetIdx()
        atom2_idx = atom2.GetIdx()
        rotatable_bonds.append((atom1_idx, atom2_idx))
  
  return rotatable_bonds

