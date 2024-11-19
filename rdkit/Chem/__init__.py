#
#  Copyright (C) 2000-2017  greg Landrum and other RDKit contributors
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" A module for molecules and stuff

 see Chem/index.html in the doc tree for documentation

"""
from rdkit import DataStructs, RDConfig, rdBase
from rdkit.Chem import rdchem
from rdkit.Geometry import rdGeometry

_HasSubstructMatchStr = rdchem._HasSubstructMatchStr
from rdkit.Chem.inchi import *
from rdkit.Chem.rdchem import *
from rdkit.Chem.rdCIPLabeler import *
from rdkit.Chem.rdmolfiles import *
from rdkit.Chem.rdmolops import *

try:
  # This is an optional component of the build
  from rdkit.Chem.rdMolInterchange import *
except ImportError:
  pass

# Coordgen needs to know where its template file is.
# The default install puts it in RDDataDir
try:
  from rdkit.Chem import rdCoordGen
except ImportError:
  pass
else:
  templDir = RDConfig.RDDataDir
  if templDir[-1] != '/':
    templDir += '/'
  rdCoordGen.SetDefaultTemplateFileDir(templDir)


# Github Issue #6208: boost::python iterators are slower than they should
# (cause is that boost::python throws exceptions as actual C++ exceptions)
class _GetRDKitObjIterator:

  def _sizeCalc(self):
    raise NotImplementedError()

  def _getRDKitItem(self, i):
    raise NotImplementedError()

  def __init__(self, mol):
    self._mol = mol
    self._pos = 0
    self._size = self._sizeCalc()

  def __len__(self):
    if self._sizeCalc() != self._size:
      raise RuntimeError('size changed during iteration')
    return self._size

  def __getitem__(self, i):
    if i < 0 or i >= len(self):
      raise IndexError('index out of range')
    return self._getRDKitItem(i)

  def __next__(self):
    if self._pos >= len(self):
      raise StopIteration

    ret = self[self._pos]
    self._pos += 1
    return ret

  def __iter__(self):
    for i in range(0, len(self)):
      self._pos = i
      yield self[i]
    self._pos = self._size


class _GetAtomsIterator(_GetRDKitObjIterator):

  def _sizeCalc(self):
    return self._mol.GetNumAtoms()

  def _getRDKitItem(self, i):
    return self._mol.GetAtomWithIdx(i)


class _GetBondsIterator(_GetRDKitObjIterator):

  def _sizeCalc(self):
    return self._mol.GetNumBonds()

  def _getRDKitItem(self, i):
    return self._mol.GetBondWithIdx(i)


rdchem.Mol.GetAtoms = lambda self: _GetAtomsIterator(self)
rdchem.Mol.GetAtoms.__doc__ = """returns an iterator over the atoms in the molecule"""
rdchem.Mol.GetBonds = lambda self: _GetBondsIterator(self)
rdchem.Mol.GetBonds.__doc__ = """returns an iterator over the bonds in the molecule"""


def QuickSmartsMatch(smi, sma, unique=True, display=False):
  ''' A convenience function for quickly matching a SMARTS against a SMILES

  Arguments:
    - smi: the SMILES to match
    - sma: the SMARTS to match
    - unique: (optional) determines whether or not only unique matches are returned
    - display: (optional) IGNORED

  Returns:
    a list of list of the indices of the atoms in the molecule that match the SMARTS  
  
  '''
  m = MolFromSmiles(smi)
  p = MolFromSmarts(sma)
  res = m.GetSubstructMatches(p, unique)
  return res


def CanonSmiles(smi, useChiral=1):
  ''' A convenience function for canonicalizing SMILES

  Arguments:
    - smi: the SMILES to canonicalize
    - useChiral: (optional) determines whether or not chiral information is included in the canonicalization and SMILES

  Returns:
    the canonical SMILES

  '''
  m = MolFromSmiles(smi)
  return MolToSmiles(m, useChiral)


def SupplierFromFilename(fileN, delim='', **kwargs):
  ''' A convenience function for creating a molecule supplier from a filename 
  
  Arguments:
    - fileN: the name of the file to read from
    - delim: (optional) the delimiter to use for reading the file (only for csv and txt files)
    - kwargs: additional keyword arguments to be passed to the supplier constructor

  Returns:
    a molecule supplier

  '''
  ext = fileN.split('.')[-1].lower()
  if ext == 'sdf':
    suppl = SDMolSupplier(fileN, **kwargs)
  elif ext == 'csv':
    if not delim:
      delim = ','
    suppl = SmilesMolSupplier(fileN, delimiter=delim, **kwargs)
  elif ext == 'txt':
    if not delim:
      delim = '\t'
    suppl = SmilesMolSupplier(fileN, delimiter=delim, **kwargs)
  elif ext == 'tdt':
    suppl = TDTMolSupplier(fileN, delimiter=delim, **kwargs)
  else:
    raise ValueError("unrecognized extension: %s" % ext)

  return suppl


def FindMolChiralCenters(mol, force=True, includeUnassigned=False, includeCIP=True,
                         useLegacyImplementation=None):
  """ returns information about the chiral centers in a molecule

  Arguments:
    - mol: the molecule to work with
    - force: (optional) if True, stereochemistry will be assigned even if it has been already
    - includeUnassigned: (optional) if True, unassigned stereo centers will be included in the output
    - includeCIP: (optional) if True, the CIP code for each chiral center will be included in the output
    - useLegacyImplementation: (optional) if True, the legacy stereochemistry perception code will be used

  Returns:
    a list of tuples of the form (atomId, CIPCode)

    >>> from rdkit import Chem
    >>> mol = Chem.MolFromSmiles('[C@H](Cl)(F)Br')
    >>> Chem.FindMolChiralCenters(mol)
    [(0, 'R')]
    >>> mol = Chem.MolFromSmiles('[C@@H](Cl)(F)Br')
    >>> Chem.FindMolChiralCenters(mol)
    [(0, 'S')]

    >>> Chem.FindMolChiralCenters(Chem.MolFromSmiles('CCC'))
    []

    By default unassigned stereo centers are not reported:

    >>> mol = Chem.MolFromSmiles('C[C@H](F)C(F)(Cl)Br')
    >>> Chem.FindMolChiralCenters(mol,force=True)
    [(1, 'S')]

    but this can be changed:

    >>> Chem.FindMolChiralCenters(mol,force=True,includeUnassigned=True)
    [(1, 'S'), (3, '?')]

    The handling of unassigned stereocenters for dependent stereochemistry is not correct 
    using the legacy implementation:

    >>> Chem.FindMolChiralCenters(Chem.MolFromSmiles('C1CC(C)C(C)C(C)C1'),includeUnassigned=True)
    [(2, '?'), (6, '?')]
    >>> Chem.FindMolChiralCenters(Chem.MolFromSmiles('C1C[C@H](C)C(C)[C@H](C)C1'),includeUnassigned=True)
    [(2, 'S'), (4, '?'), (6, 'R')]

    But works with the new implementation:

    >>> Chem.FindMolChiralCenters(Chem.MolFromSmiles('C1CC(C)C(C)C(C)C1'),includeUnassigned=True, useLegacyImplementation=False)
    [(2, '?'), (4, '?'), (6, '?')]

    Note that the new implementation also gets the correct descriptors for para-stereochemistry:

    >>> Chem.FindMolChiralCenters(Chem.MolFromSmiles('C1C[C@H](C)[C@H](C)[C@H](C)C1'),useLegacyImplementation=False)
    [(2, 'S'), (4, 's'), (6, 'R')]

    With the new implementation, if you don't care about the CIP labels of stereocenters, you can save
    some time by disabling those:

    >>> Chem.FindMolChiralCenters(Chem.MolFromSmiles('C1C[C@H](C)[C@H](C)[C@H](C)C1'), includeCIP=False, useLegacyImplementation=False)
    [(2, 'Tet_CCW'), (4, 'Tet_CCW'), (6, 'Tet_CCW')]

  """
  origUseLegacyVal = GetUseLegacyStereoPerception()
  if useLegacyImplementation is None:
    useLegacyImplementation = origUseLegacyVal
  SetUseLegacyStereoPerception(useLegacyImplementation)
  try:
    if useLegacyImplementation:
      AssignStereochemistry(mol, force=force, flagPossibleStereoCenters=includeUnassigned)
      centers = []
      for atom in mol.GetAtoms():
        if atom.HasProp('_CIPCode'):
          centers.append((atom.GetIdx(), atom.GetProp('_CIPCode')))
        elif includeUnassigned and atom.HasProp('_ChiralityPossible'):
          centers.append((atom.GetIdx(), '?'))
    else:
      centers = []
      itms = FindPotentialStereo(mol)
      if includeCIP:
        atomsToLabel = []
        bondsToLabel = []
        for si in itms:
          if si.type == StereoType.Atom_Tetrahedral:
            atomsToLabel.append(si.centeredOn)
          elif si.type == StereoType.Bond_Double:
            bondsToLabel.append(si.centeredOn)
        AssignCIPLabels(mol, atomsToLabel=atomsToLabel, bondsToLabel=bondsToLabel)
      for si in itms:
        if si.type == StereoType.Atom_Tetrahedral and (includeUnassigned or si.specified
                                                       == StereoSpecified.Specified):
          idx = si.centeredOn
          atm = mol.GetAtomWithIdx(idx)
          if includeCIP and atm.HasProp("_CIPCode"):
            code = atm.GetProp("_CIPCode")
          else:
            if si.specified:
              code = str(si.descriptor)
            else:
              code = '?'
              atm.SetIntProp('_ChiralityPossible', 1)
          centers.append((idx, code))
  finally:
    SetUseLegacyStereoPerception(origUseLegacyVal)
  return centers
