# $Id$
#
# Copyright (C) 2003-2006 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" Calculation of topological/topochemical descriptors.



"""

import math

import numpy

from rdkit import Chem
from rdkit.Chem import Graphs, rdchem, rdMolDescriptors
from rdkit.ML.InfoTheory import entropy

ptable = Chem.GetPeriodicTable()

_log2val = math.log(2)


def _VertexDegrees(mat, onlyOnes=0):
  """  *Internal Use Only*

  this is just a row sum of the matrix... simple, neh?

  """
  if not onlyOnes:
    res = sum(mat)
  else:
    res = sum(numpy.equal(mat, 1))
  return res


def _NumAdjacencies(mol, dMat):
  """  *Internal Use Only*

  """
  res = mol.GetNumBonds()
  return res


def _GetCountDict(arr):
  """  *Internal Use Only*

  """
  res = {}
  for v in arr:
    res[v] = res.get(v, 0) + 1
  return res


# WARNING: this data should probably go somewhere else...
hallKierAlphas = {
  'Br': [None, None, 0.48],
  'C': [-0.22, -0.13, 0.0],
  'Cl': [None, None, 0.29],
  'F': [None, None, -0.07],
  'H': [0.0, 0.0, 0.0],
  'I': [None, None, 0.73],
  'N': [-0.29, -0.2, -0.04],
  'O': [None, -0.2, -0.04],
  'P': [None, 0.3, 0.43],
  'S': [None, 0.22, 0.35]
}


def _pyHallKierAlpha(m):
  """ calculate the Hall-Kier alpha value for a molecule

   From equations (58) of Rev. Comp. Chem. vol 2, 367-422, (1991)

  """
  alphaSum = 0.0
  rC = ptable.GetRb0(6)
  for atom in m.GetAtoms():
    atNum = atom.GetAtomicNum()
    if not atNum:
      continue
    symb = atom.GetSymbol()
    alphaV = hallKierAlphas.get(symb, None)
    if alphaV is not None:
      hyb = atom.GetHybridization() - 2
      if (hyb < len(alphaV)):
        alpha = alphaV[hyb]
        if alpha is None:
          alpha = alphaV[-1]
      else:
        alpha = alphaV[-1]
    else:
      rA = ptable.GetRb0(atNum)
      alpha = rA / rC - 1
    # print(atom.GetIdx(), atom.GetSymbol(), alpha)
    alphaSum += alpha
  return alphaSum


# HallKierAlpha.version="1.0.2"


def Ipc(mol, avg=False, dMat=None, forceDMat=False):
  """This returns the information content of the coefficients of the characteristic
    polynomial of the adjacency matrix of a hydrogen-suppressed graph of a molecule.

    'avg = True' returns the information content divided by the total population.

    From Eq 6 of D. Bonchev & N. Trinajstic, J. Chem. Phys. vol 67, 4517-4533 (1977)

  """
  if forceDMat or dMat is None:
    if forceDMat:
      dMat = Chem.GetDistanceMatrix(mol, 0)
      mol._adjMat = dMat
    else:
      try:
        dMat = mol._adjMat
      except AttributeError:
        dMat = Chem.GetDistanceMatrix(mol, 0)
        mol._adjMat = dMat

  adjMat = numpy.equal(dMat, 1)
  cPoly = abs(Graphs.CharacteristicPolynomial(mol, adjMat))
  if avg:
    return entropy.InfoEntropy(cPoly)
  else:
    return sum(cPoly) * entropy.InfoEntropy(cPoly)


Ipc.version = "1.0.0"


def AvgIpc(mol, dMat=None, forceDMat=False):
  """This returns the average information content of the coefficients of the characteristic
    polynomial of the adjacency matrix of a hydrogen-suppressed graph of a molecule.

    From Eq 7 of D. Bonchev & N. Trinajstic, J. Chem. Phys. vol 67, 4517-4533 (1977)

  """
  return Ipc(mol, avg=True, dMat=dMat, forceDMat=forceDMat)


AvgIpc.version = "1.0.0"


def _pyKappa1(mol):
  """ Hall-Kier Kappa1 value

   From equations (58) and (59) of Rev. Comp. Chem. vol 2, 367-422, (1991)

  """
  P1 = mol.GetNumBonds(1)
  A = mol.GetNumHeavyAtoms()
  alpha = HallKierAlpha(mol)
  denom = P1 + alpha
  if denom:
    kappa = (A + alpha) * (A + alpha - 1)**2 / denom**2
  else:
    kappa = 0.0
  return kappa


# Kappa1.version="1.0.0"


def _pyKappa2(mol):
  """  Hall-Kier Kappa2 value

   From equations (58) and (60) of Rev. Comp. Chem. vol 2, 367-422, (1991)

  """
  P2 = len(Chem.FindAllPathsOfLengthN(mol, 2))
  A = mol.GetNumHeavyAtoms()
  alpha = HallKierAlpha(mol)
  denom = (P2 + alpha)**2
  if denom:
    kappa = (A + alpha - 1) * (A + alpha - 2)**2 / denom
  else:
    kappa = 0
  return kappa


# Kappa2.version="1.0.0"


def _pyKappa3(mol):
  """  Hall-Kier Kappa3 value

   From equations (58), (61) and (62) of Rev. Comp. Chem. vol 2, 367-422, (1991)

  """
  P3 = len(Chem.FindAllPathsOfLengthN(mol, 3))
  A = mol.GetNumHeavyAtoms()
  alpha = HallKierAlpha(mol)
  denom = (P3 + alpha)**2
  if denom:
    if A % 2 == 1:
      kappa = (A + alpha - 1) * (A + alpha - 3)**2 / denom
    else:
      kappa = (A + alpha - 2) * (A + alpha - 3)**2 / denom
  else:
    kappa = 0
  return kappa


# Kappa3.version="1.0.0"

HallKierAlpha = lambda x: rdMolDescriptors.CalcHallKierAlpha(x)
HallKierAlpha.version = rdMolDescriptors._CalcHallKierAlpha_version
Kappa1 = lambda x: rdMolDescriptors.CalcKappa1(x)
Kappa1.version = rdMolDescriptors._CalcKappa1_version
Kappa2 = lambda x: rdMolDescriptors.CalcKappa2(x)
Kappa2.version = rdMolDescriptors._CalcKappa2_version
Kappa3 = lambda x: rdMolDescriptors.CalcKappa3(x)
Kappa3.version = rdMolDescriptors._CalcKappa3_version


def Chi0(mol):
  """ From equations (1),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, (1991)

  """
  deltas = [x.GetDegree() for x in mol.GetAtoms()]
  while 0 in deltas:
    deltas.remove(0)
  deltas = numpy.array(deltas, 'd')
  res = sum(numpy.sqrt(1. / deltas))
  return res


Chi0.version = "1.0.0"


def Chi1(mol):
  """ From equations (1),(11) and (12) of Rev. Comp. Chem. vol 2, 367-422, (1991)

  """
  c1s = [x.GetBeginAtom().GetDegree() * x.GetEndAtom().GetDegree() for x in mol.GetBonds()]
  while 0 in c1s:
    c1s.remove(0)
  c1s = numpy.array(c1s, 'd')
  res = sum(numpy.sqrt(1. / c1s))
  return res


Chi1.version = "1.0.0"


def _nVal(atom):
  return ptable.GetNOuterElecs(atom.GetAtomicNum()) - atom.GetTotalNumHs()


def _hkDeltas(mol, skipHs=1):
  global ptable
  res = []
  if hasattr(mol, '_hkDeltas') and mol._hkDeltas is not None:
    return mol._hkDeltas
  for atom in mol.GetAtoms():
    n = atom.GetAtomicNum()
    if n > 1:
      nV = ptable.GetNOuterElecs(n)
      nHs = atom.GetTotalNumHs()
      if n <= 10:
        # first row
        res.append(float(nV - nHs))
      else:
        # second row and up
        res.append(float(nV - nHs) / float(n - nV - 1))
    elif n == 1:
      if not skipHs:
        res.append(0.0)
    else:
      res.append(0.0)
  mol._hkDeltas = res
  return res


def _pyChi0v(mol):
  """  From equations (5),(9) and (10) of Rev. Comp. Chem. vol 2, 367-422, (1991)

  """
  deltas = _hkDeltas(mol)
  while 0 in deltas:
    deltas.remove(0)
  mol._hkDeltas = None
  res = sum(numpy.sqrt(1. / numpy.array(deltas)))
  return res


def _pyChi1v(mol):
  """  From equations (5),(11) and (12) of Rev. Comp. Chem. vol 2, 367-422, (1991)

  """
  deltas = numpy.array(_hkDeltas(mol, skipHs=0))
  res = 0.0
  for bond in mol.GetBonds():
    v = deltas[bond.GetBeginAtomIdx()] * deltas[bond.GetEndAtomIdx()]
    if v != 0.0:
      res += numpy.sqrt(1. / v)
  return res


def _pyChiNv_(mol, order=2):
  """  From equations (5),(15) and (16) of Rev. Comp. Chem. vol 2, 367-422, (1991)

  **NOTE**: because the current path finding code does, by design,
  detect rings as paths (e.g. in C1CC1 there is *1* atom path of
  length 3), values of ChiNv with N >= 3 may give results that differ
  from those provided by the old code in molecules that have rings of
  size 3.

  """
  deltas = numpy.array([(1. / numpy.sqrt(hkd) if hkd != 0.0 else 0.0)
                        for hkd in _hkDeltas(mol, skipHs=0)])
  accum = 0.0
  for path in Chem.FindAllPathsOfLengthN(mol, order + 1, useBonds=0):
    accum += numpy.prod(deltas[numpy.array(path)])
  return accum


def _pyChi2v(mol):
  """ From equations (5),(15) and (16) of Rev. Comp. Chem. vol 2, 367-422, (1991)

  """
  return _pyChiNv_(mol, 2)


def _pyChi3v(mol):
  """ From equations (5),(15) and (16) of Rev. Comp. Chem. vol 2, 367-422, (1991)

  """
  return _pyChiNv_(mol, 3)


def _pyChi4v(mol):
  """ From equations (5),(15) and (16) of Rev. Comp. Chem. vol 2, 367-422, (1991)

  **NOTE**: because the current path finding code does, by design,
  detect rings as paths (e.g. in C1CC1 there is *1* atom path of
  length 3), values of Chi4v may give results that differ from those
  provided by the old code in molecules that have 3 rings.

  """
  return _pyChiNv_(mol, 4)


def _pyChi0n(mol):
  """  Similar to Hall Kier Chi0v, but uses nVal instead of valence
  This makes a big difference after we get out of the first row.

  """
  deltas = [_nVal(x) for x in mol.GetAtoms()]
  while deltas.count(0):
    deltas.remove(0)
  deltas = numpy.array(deltas, 'd')
  res = sum(numpy.sqrt(1. / deltas))
  return res


def _pyChi1n(mol):
  """  Similar to Hall Kier Chi1v, but uses nVal instead of valence

  """
  delts = numpy.array([_nVal(x) for x in mol.GetAtoms()], 'd')
  res = 0.0
  for bond in mol.GetBonds():
    v = delts[bond.GetBeginAtomIdx()] * delts[bond.GetEndAtomIdx()]
    if v != 0.0:
      res += numpy.sqrt(1. / v)
  return res


def _pyChiNn_(mol, order=2):
  """  Similar to Hall Kier ChiNv, but uses nVal instead of valence
  This makes a big difference after we get out of the first row.

  **NOTE**: because the current path finding code does, by design,
  detect rings as paths (e.g. in C1CC1 there is *1* atom path of
  length 3), values of ChiNn with N >= 3 may give results that differ
  from those provided by the old code in molecules that have rings of
  size 3.

  """
  nval = [_nVal(x) for x in mol.GetAtoms()]
  deltas = numpy.array([(1. / numpy.sqrt(x) if x else 0.0) for x in nval])
  accum = 0.0
  for path in Chem.FindAllPathsOfLengthN(mol, order + 1, useBonds=0):
    accum += numpy.prod(deltas[numpy.array(path)])
  return accum


def _pyChi2n(mol):
  """  Similar to Hall Kier Chi2v, but uses nVal instead of valence
  This makes a big difference after we get out of the first row.

  """
  return _pyChiNn_(mol, 2)


def _pyChi3n(mol):
  """  Similar to Hall Kier Chi3v, but uses nVal instead of valence
  This makes a big difference after we get out of the first row.

  """
  return _pyChiNn_(mol, 3)


def _pyChi4n(mol):
  """  Similar to Hall Kier Chi4v, but uses nVal instead of valence
  This makes a big difference after we get out of the first row.


  **NOTE**: because the current path finding code does, by design,
  detect rings as paths (e.g. in C1CC1 there is *1* atom path of
  length 3), values of Chi4n may give results that differ from those
  provided by the old code in molecules that have 3 rings.

  """
  return _pyChiNn_(mol, 4)


Chi0v = lambda x: rdMolDescriptors.CalcChi0v(x)
Chi0v.version = rdMolDescriptors._CalcChi0v_version
Chi1v = lambda x: rdMolDescriptors.CalcChi1v(x)
Chi1v.version = rdMolDescriptors._CalcChi1v_version
Chi2v = lambda x: rdMolDescriptors.CalcChi2v(x)
Chi2v.version = rdMolDescriptors._CalcChi2v_version
Chi3v = lambda x: rdMolDescriptors.CalcChi3v(x)
Chi3v.version = rdMolDescriptors._CalcChi3v_version
Chi4v = lambda x: rdMolDescriptors.CalcChi4v(x)
Chi4v.version = rdMolDescriptors._CalcChi4v_version
ChiNv_ = lambda x, y: rdMolDescriptors.CalcChiNv(x, y)
ChiNv_.version = rdMolDescriptors._CalcChiNv_version

Chi0n = lambda x: rdMolDescriptors.CalcChi0n(x)
Chi0n.version = rdMolDescriptors._CalcChi0n_version
Chi1n = lambda x: rdMolDescriptors.CalcChi1n(x)
Chi1n.version = rdMolDescriptors._CalcChi1n_version
Chi2n = lambda x: rdMolDescriptors.CalcChi2n(x)
Chi2n.version = rdMolDescriptors._CalcChi2n_version
Chi3n = lambda x: rdMolDescriptors.CalcChi3n(x)
Chi3n.version = rdMolDescriptors._CalcChi3n_version
Chi4n = lambda x: rdMolDescriptors.CalcChi4n(x)
Chi4n.version = rdMolDescriptors._CalcChi4n_version
ChiNn_ = lambda x, y: rdMolDescriptors.CalcChiNn(x, y)
ChiNn_.version = rdMolDescriptors._CalcChiNn_version


def BalabanJ(mol, dMat=None, forceDMat=0):
  """ Calculate Balaban's J value for a molecule

  **Arguments**

    - mol: a molecule

    - dMat: (optional) a distance/adjacency matrix for the molecule, if this
      is not provide, one will be calculated

    - forceDMat: (optional) if this is set, the distance/adjacency matrix
      will be recalculated regardless of whether or not _dMat_ is provided
      or the molecule already has one

  **Returns**

    - a float containing the J value

  We follow the notation of Balaban's paper:
    Chem. Phys. Lett. vol 89, 399-404, (1982)

  """
  # if no dMat is passed in, calculate one ourselves
  if forceDMat or dMat is None:
    if forceDMat:
      # FIX: should we be using atom weights here or not?
      dMat = Chem.GetDistanceMatrix(mol, useBO=1, useAtomWts=0, force=1)
      mol._balabanMat = dMat
      adjMat = Chem.GetAdjacencyMatrix(mol, useBO=0, emptyVal=0, force=0, prefix="NoBO")
      mol._adjMat = adjMat
    else:
      try:
        # first check if the molecule already has one
        dMat = mol._balabanMat
      except AttributeError:
        # nope, gotta calculate one
        dMat = Chem.GetDistanceMatrix(mol, useBO=1, useAtomWts=0, force=0, prefix="Balaban")
        # now store it
        mol._balabanMat = dMat
      try:
        adjMat = mol._adjMat
      except AttributeError:
        adjMat = Chem.GetAdjacencyMatrix(mol, useBO=0, emptyVal=0, force=0, prefix="NoBO")
        mol._adjMat = adjMat
  else:
    adjMat = Chem.GetAdjacencyMatrix(mol, useBO=0, emptyVal=0, force=0, prefix="NoBO")

  s = _VertexDegrees(dMat)
  q = _NumAdjacencies(mol, dMat)
  n = mol.GetNumAtoms()
  mu = q - n + 1

  sum_ = 0.
  nS = len(s)
  for i in range(nS):
    si = s[i]
    for j in range(i, nS):
      if adjMat[i, j] == 1:
        sum_ += 1. / numpy.sqrt(si * s[j])

  if mu + 1 != 0:
    J = float(q) / float(mu + 1) * sum_
  else:
    J = 0

  return J


BalabanJ.version = "1.0.0"

# ------------------------------------------------------------------------
#
# Start block of BertzCT stuff.
#


def _AssignSymmetryClasses(mol, vdList, bdMat, forceBDMat, numAtoms, cutoff):
  """
     Used by BertzCT

     vdList: the number of neighbors each atom has
     bdMat: "balaban" distance matrix

  """
  if forceBDMat:
    bdMat = Chem.GetDistanceMatrix(mol, useBO=1, useAtomWts=0, force=1, prefix="Balaban")
    mol._balabanMat = bdMat

  keysSeen = []
  symList = [0] * numAtoms
  for i in range(numAtoms):
    tmpList = bdMat[i].tolist()
    tmpList.sort()
    theKey = tuple(['%.4f' % x for x in tmpList[:cutoff]])
    try:
      idx = keysSeen.index(theKey)
    except ValueError:
      idx = len(keysSeen)
      keysSeen.append(theKey)
    symList[i] = idx + 1
  return tuple(symList)


def _LookUpBondOrder(atom1Id, atom2Id, bondDic):
  """
     Used by BertzCT
  """
  if atom1Id < atom2Id:
    theKey = (atom1Id, atom2Id)
  else:
    theKey = (atom2Id, atom1Id)
  tmp = bondDic[theKey]
  if tmp == Chem.BondType.AROMATIC:
    tmp = 1.5
  else:
    tmp = float(tmp)
  # tmp = int(tmp)
  return tmp


def _CalculateEntropies(connectionDict, atomTypeDict, numAtoms):
  """
     Used by BertzCT
  """
  connectionList = list(connectionDict.values())
  totConnections = sum(connectionList)
  connectionIE = totConnections * (entropy.InfoEntropy(numpy.array(connectionList)) +
                                   math.log(totConnections) / _log2val)
  atomTypeList = list(atomTypeDict.values())
  atomTypeIE = numAtoms * entropy.InfoEntropy(numpy.array(atomTypeList))
  return atomTypeIE + connectionIE


def _CreateBondDictEtc(mol, numAtoms):
  """ _Internal Use Only_
     Used by BertzCT

  """
  bondDict = {}
  nList = [None] * numAtoms
  vdList = [0] * numAtoms
  for aBond in mol.GetBonds():
    atom1 = aBond.GetBeginAtomIdx()
    atom2 = aBond.GetEndAtomIdx()
    if atom1 > atom2:
      atom2, atom1 = atom1, atom2
    if not aBond.GetIsAromatic():
      bondDict[(atom1, atom2)] = aBond.GetBondType()
    else:
      # mark Kekulized systems as aromatic
      bondDict[(atom1, atom2)] = Chem.BondType.AROMATIC
    if nList[atom1] is None:
      nList[atom1] = [atom2]
    elif atom2 not in nList[atom1]:
      nList[atom1].append(atom2)
    if nList[atom2] is None:
      nList[atom2] = [atom1]
    elif atom1 not in nList[atom2]:
      nList[atom2].append(atom1)

  for i, element in enumerate(nList):
    try:
      element.sort()
      vdList[i] = len(element)
    except Exception:
      vdList[i] = 0
  return bondDict, nList, vdList


def BertzCT(mol, cutoff=100, dMat=None, forceDMat=1):
  """ A topological index meant to quantify "complexity" of molecules.

     Consists of a sum of two terms, one representing the complexity
     of the bonding, the other representing the complexity of the
     distribution of heteroatoms.

     From S. H. Bertz, J. Am. Chem. Soc., vol 103, 3599-3601 (1981)

     "cutoff" is an integer value used to limit the computational
     expense.  A cutoff value tells the program to consider vertices
     topologically identical if their distance vectors (sets of
     distances to all other vertices) are equal out to the "cutoff"th
     nearest-neighbor.

     **NOTE**  The original implementation had the following comment:
         > this implementation treats aromatic rings as the
         > corresponding Kekule structure with alternating bonds,
         > for purposes of counting "connections".
       Upon further thought, this is the WRONG thing to do.  It
        results in the possibility of a molecule giving two different
        CT values depending on the kekulization.  For example, in the
        old implementation, these two SMILES:
           CC2=CN=C1C3=C(C(C)=C(C=N3)C)C=CC1=C2C
           CC3=CN=C2C1=NC=C(C)C(C)=C1C=CC2=C3C
        which correspond to differentk kekule forms, yield different
        values.
       The new implementation uses consistent (aromatic) bond orders
        for aromatic bonds.

       THIS MEANS THAT THIS IMPLEMENTATION IS NOT BACKWARDS COMPATIBLE.

       Any molecule containing aromatic rings will yield different
       values with this implementation.  The new behavior is the correct
       one, so we're going to live with the breakage.

     **NOTE** this barfs if the molecule contains a second (or
       nth) fragment that is one atom.

  """
  atomTypeDict = {}
  connectionDict = {}
  numAtoms = mol.GetNumAtoms()
  if forceDMat or dMat is None:
    if forceDMat:
      # nope, gotta calculate one
      dMat = Chem.GetDistanceMatrix(mol, useBO=0, useAtomWts=0, force=1)
      mol._adjMat = dMat
    else:
      try:
        dMat = mol._adjMat
      except AttributeError:
        dMat = Chem.GetDistanceMatrix(mol, useBO=0, useAtomWts=0, force=1)
        mol._adjMat = dMat

  if numAtoms < 2:
    return 0

  bondDict, neighborList, vdList = _CreateBondDictEtc(mol, numAtoms)
  symmetryClasses = _AssignSymmetryClasses(mol, vdList, dMat, forceDMat, numAtoms, cutoff)
  # print('Symmm Classes:',symmetryClasses)
  for atomIdx in range(numAtoms):
    hingeAtomNumber = mol.GetAtomWithIdx(atomIdx).GetAtomicNum()
    atomTypeDict[hingeAtomNumber] = atomTypeDict.get(hingeAtomNumber, 0) + 1

    hingeAtomClass = symmetryClasses[atomIdx]
    numNeighbors = vdList[atomIdx]
    for i in range(numNeighbors):
      neighbor_iIdx = neighborList[atomIdx][i]
      NiClass = symmetryClasses[neighbor_iIdx]
      bond_i_order = _LookUpBondOrder(atomIdx, neighbor_iIdx, bondDict)
      # print('\t',atomIdx,i,hingeAtomClass,NiClass,bond_i_order)
      if (bond_i_order > 1) and (neighbor_iIdx > atomIdx):
        numConnections = bond_i_order * (bond_i_order - 1) / 2
        connectionKey = (min(hingeAtomClass, NiClass), max(hingeAtomClass, NiClass))
        connectionDict[connectionKey] = connectionDict.get(connectionKey, 0) + numConnections

      for j in range(i + 1, numNeighbors):
        neighbor_jIdx = neighborList[atomIdx][j]
        NjClass = symmetryClasses[neighbor_jIdx]
        bond_j_order = _LookUpBondOrder(atomIdx, neighbor_jIdx, bondDict)
        numConnections = bond_i_order * bond_j_order
        connectionKey = (min(NiClass, NjClass), hingeAtomClass, max(NiClass, NjClass))
        connectionDict[connectionKey] = connectionDict.get(connectionKey, 0) + numConnections

  if not connectionDict:
    connectionDict = {'a': 1}

  return _CalculateEntropies(connectionDict, atomTypeDict, numAtoms)


BertzCT.version = "2.0.0"
# Recent Revisions:
#  1.0.0 -> 2.0.0:
#    - force distance matrix updates properly (Fixed as part of Issue 125)
#    - handle single-atom fragments (Issue 136)

#
# End block of BertzCT stuff.
#
# ------------------------------------------------------------------------
