#  Copyright (C) 2014-2024 Sereina Riniker and other RDKit contributors
#
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" Torsion Fingerprints (Deviation) (TFD)
    According to a paper from Schulz-Gasch et al., JCIM, 52, 1499-1512 (2012).

"""
import math
import os

from rdkit import Chem, Geometry, RDConfig, rdBase
from rdkit.Chem import rdchem, rdMolDescriptors
from rdkit.Chem import rdFingerprintGenerator

def _doMatch(inv, atoms):
  """ Helper function to check if all atoms in the list are the same
      
      Arguments:
      - inv:    atom invariants (used to define equivalence of atoms)
      - atoms:  list of atoms to check

      Return: boolean
  """
  match = True
  for i in range(len(atoms) - 1):
    for j in range(i + 1, len(atoms)):
      if (inv[atoms[i].GetIdx()] != inv[atoms[j].GetIdx()]):
        match = False
        return match
  return match


def _doNotMatch(inv, atoms):
  """ Helper function to check if all atoms in the list are NOT the same
      
      Arguments:
      - inv:    atom invariants (used to define equivalence of atoms)
      - atoms:  list of atoms to check

      Return: boolean
  """
  match = True
  for i in range(len(atoms) - 1):
    for j in range(i + 1, len(atoms)):
      if (inv[atoms[i].GetIdx()] == inv[atoms[j].GetIdx()]):
        match = False
        return match
  return match


def _doMatchExcept1(inv, atoms):
  """ Helper function to check if two atoms in the list are the same, 
      and one not
      Note: Works only for three atoms
      
      Arguments:
      - inv:    atom invariants (used to define equivalence of atoms)
      - atoms:  list of atoms to check

      Return: atom that is different
  """
  if len(atoms) != 3:
    raise ValueError("Number of atoms must be three")
  a1 = atoms[0].GetIdx()
  a2 = atoms[1].GetIdx()
  a3 = atoms[2].GetIdx()
  if (inv[a1] == inv[a2] and inv[a1] != inv[a3] and inv[a2] != inv[a3]):
    return atoms[2]
  elif (inv[a1] != inv[a2] and inv[a1] == inv[a3] and inv[a2] != inv[a3]):
    return atoms[1]
  elif (inv[a1] != inv[a2] and inv[a1] != inv[a3] and inv[a2] == inv[a3]):
    return atoms[0]
  return None


def _getAtomInvariantsWithRadius(mol, radius):
  """ Helper function to calculate the atom invariants for each atom 
      with a given radius

      Arguments:
      - mol:    the molecule of interest
      - radius: the radius for the Morgan fingerprint

      Return: list of atom invariants
  """
  inv = []
  fpg = rdFingerprintGenerator.GetMorganGenerator(radius,includeRedundantEnvironments=True)
  ao = rdFingerprintGenerator.AdditionalOutput()
  ao.AllocateBitInfoMap()

  fp = fpg.GetSparseCountFingerprint(mol, additionalOutput=ao)
  bim = ao.GetBitInfoMap()
  inv = [0]*mol.GetNumAtoms()
  for k,vs in bim.items():
    for (aid,r) in vs:
      if r == radius:
        inv[aid] = k
  return inv


def _getHeavyAtomNeighbors(atom1, aid2=-1):
  """ Helper function to calculate the number of heavy atom neighbors.

      Arguments:
      - atom1:    the atom of interest
      - aid2:     atom index that should be excluded from neighbors (default: none)

      Return: a list of heavy atom neighbors of the given atom
  """
  if aid2 < 0:
    return [n for n in atom1.GetNeighbors() if n.GetSymbol() != 'H']
  else:
    return [n for n in atom1.GetNeighbors() if (n.GetSymbol() != 'H' and n.GetIdx() != aid2)]


def _getIndexforTorsion(neighbors, inv):
  """ Helper function to calculate the index of the reference atom for 
      a given atom

      Arguments:
      - neighbors:  list of the neighbors of the atom
      - inv:        atom invariants

      Return: list of atom indices as reference for torsion
  """
  if len(neighbors) == 1:  # atom has only one neighbor
    return [neighbors[0]]
  elif _doMatch(inv, neighbors):  # atom has all symmetric neighbors
    return neighbors
  elif _doNotMatch(inv, neighbors):  # atom has all different neighbors
    # sort by atom inv and simply use the first neighbor
    neighbors = sorted(neighbors, key=lambda x: inv[x.GetIdx()])
    return [neighbors[0]]
  elif len(neighbors) == 3:
    at = _doMatchExcept1(inv, neighbors)  # two neighbors the same, one different
    if at is None:
      raise ValueError("Atom neighbors are either all the same or all different")
    return [at]
  else:  # weird case
    # sort by atom inv and simply use the first neighbor
    neighbors = sorted(neighbors, key=lambda x: inv[x.GetIdx()])
    return [neighbors[0]]


def _getBondsForTorsions(mol, ignoreColinearBonds):
  """ Determine the bonds (or pair of atoms treated like a bond) for which
      torsions should be calculated.

      Arguments:
      - refmol: the molecule of interest
      - ignoreColinearBonds: if True (default), single bonds adjacent to
                             triple bonds are ignored
                             if False, alternative not-covalently bound
                             atoms are used to define the torsion
  """
  # flag the atoms that cannot be part of the centre atoms of a torsion
  # patterns: triple bonds and allenes
  patts = [Chem.MolFromSmarts(x) for x in ['*#*', '[$([C](=*)=*)]']]
  atomFlags = [0] * mol.GetNumAtoms()
  for p in patts:
    if mol.HasSubstructMatch(p):
      matches = mol.GetSubstructMatches(p)
      for match in matches:
        for a in match:
          atomFlags[a] = 1

  bonds = []
  doneBonds = [0] * mol.GetNumBonds()
  for b in mol.GetBonds():
    if b.IsInRing():
      continue
    a1 = b.GetBeginAtomIdx()
    a2 = b.GetEndAtomIdx()
    nb1 = _getHeavyAtomNeighbors(b.GetBeginAtom(), a2)
    nb2 = _getHeavyAtomNeighbors(b.GetEndAtom(), a1)
    if not doneBonds[b.GetIdx()] and (nb1 and nb2):  # no terminal bonds
      doneBonds[b.GetIdx()] = 1
      # check if atoms cannot be middle atoms
      if atomFlags[a1] or atomFlags[a2]:
        if not ignoreColinearBonds:  # search for alternative not-covalently bound atoms
          while len(nb1) == 1 and atomFlags[a1]:
            a1old = a1
            a1 = nb1[0].GetIdx()
            b = mol.GetBondBetweenAtoms(a1old, a1)
            if b.GetEndAtom().GetIdx() == a1old:
              nb1 = _getHeavyAtomNeighbors(b.GetBeginAtom(), a1old)
            else:
              nb1 = _getHeavyAtomNeighbors(b.GetEndAtom(), a1old)
            doneBonds[b.GetIdx()] = 1
          while len(nb2) == 1 and atomFlags[a2]:
            doneBonds[b.GetIdx()] = 1
            a2old = a2
            a2 = nb2[0].GetIdx()
            b = mol.GetBondBetweenAtoms(a2old, a2)
            if b.GetBeginAtom().GetIdx() == a2old:
              nb2 = _getHeavyAtomNeighbors(b.GetEndAtom(), a2old)
            else:
              nb2 = _getHeavyAtomNeighbors(b.GetBeginAtom(), a2old)
            doneBonds[b.GetIdx()] = 1
          if nb1 and nb2:
            bonds.append((a1, a2, nb1, nb2))

      else:
        bonds.append((a1, a2, nb1, nb2))
  return bonds


def CalculateTorsionLists(mol, maxDev='equal', symmRadius=2, ignoreColinearBonds=True):
  """ Calculate a list of torsions for a given molecule. For each torsion
      the four atom indices are determined and stored in a set.

      Arguments:
      - mol:      the molecule of interest
      - maxDev:   maximal deviation used for normalization
                  'equal': all torsions are normalized using 180.0 (default)
                  'spec':  each torsion is normalized using its specific
                           maximal deviation as given in the paper
      - symmRadius: radius used for calculating the atom invariants
                    (default: 2)
      - ignoreColinearBonds: if True (default), single bonds adjacent to
                             triple bonds are ignored
                             if False, alternative not-covalently bound
                             atoms are used to define the torsion

      Return: two lists of torsions: non-ring and ring torsions
  """
  if maxDev not in ['equal', 'spec']:
    raise ValueError("maxDev must be either equal or spec")
  # get non-terminal, non-cyclic bonds
  bonds = _getBondsForTorsions(mol, ignoreColinearBonds)
  # get atom invariants
  if symmRadius > 0:
    inv = _getAtomInvariantsWithRadius(mol, symmRadius)
  else:
    inv = rdMolDescriptors.GetConnectivityInvariants(mol)
  # get the torsions
  tors_list = []  # to store the atom indices of the torsions
  for a1, a2, nb1, nb2 in bonds:
    d1 = _getIndexforTorsion(nb1, inv)
    d2 = _getIndexforTorsion(nb2, inv)
    if len(d1) == 1 and len(d2) == 1:  # case 1, 2, 4, 5, 7, 10, 16, 12, 17, 19
      tors_list.append(([(d1[0].GetIdx(), a1, a2, d2[0].GetIdx())], 180.0))
    elif len(d1) == 1:  # case 3, 6, 8, 13, 20
      if len(nb2) == 2:  # two neighbors
        tors_list.append(([(d1[0].GetIdx(), a1, a2, nb.GetIdx()) for nb in d2], 90.0))
      else:  # three neighbors
        tors_list.append(([(d1[0].GetIdx(), a1, a2, nb.GetIdx()) for nb in d2], 60.0))
    elif len(d2) == 1:  # case 3, 6, 8, 13, 20
      if len(nb1) == 2:
        tors_list.append(([(nb.GetIdx(), a1, a2, d2[0].GetIdx()) for nb in d1], 90.0))
      else:  # three neighbors
        tors_list.append(([(nb.GetIdx(), a1, a2, d2[0].GetIdx()) for nb in d1], 60.0))
    else:  # both symmetric
      tmp = []
      for n1 in d1:
        for n2 in d2:
          tmp.append((n1.GetIdx(), a1, a2, n2.GetIdx()))
      if len(nb1) == 2 and len(nb2) == 2:  # case 9
        tors_list.append((tmp, 90.0))
      elif len(nb1) == 3 and len(nb2) == 3:  # case 21
        tors_list.append((tmp, 60.0))
      else:  # case 15
        tors_list.append((tmp, 30.0))
  # maximal possible deviation for non-cyclic bonds
  if maxDev == 'equal':
    tors_list = [(t, 180.0) for t, d in tors_list]
  # rings
  rings = Chem.GetSymmSSSR(mol)
  tors_list_rings = []
  for r in rings:
    # get the torsions
    tmp = []
    num = len(r)
    if 14 <= num:
      maxdev = 180.0
    else:
      maxdev = 180.0 * math.exp(-0.025 * (num - 14) * (num - 14))
    for i in range(len(r)):
      tmp.append((r[i], r[(i + 1) % num], r[(i + 2) % num], r[(i + 3) % num]))
    tors_list_rings.append((tmp, maxdev))
  return tors_list, tors_list_rings


def _getTorsionAtomPositions(atoms, conf):
  """ Helper function to retrieve the coordinates of the four atoms
      in a torsion

      Arguments:
      - atoms:   list with the four atoms
      - conf:    conformation of the molecule

      Return: Point3D objects of the four atoms
  """
  if len(atoms) != 4:
    raise ValueError("List must contain exactly four atoms")
  p1 = conf.GetAtomPosition(atoms[0])
  p2 = conf.GetAtomPosition(atoms[1])
  p3 = conf.GetAtomPosition(atoms[2])
  p4 = conf.GetAtomPosition(atoms[3])
  return p1, p2, p3, p4


def CalculateTorsionAngles(mol, tors_list, tors_list_rings, confId=-1):
  """ Calculate the torsion angles for a list of non-ring and 
      a list of ring torsions.

      Arguments:
      - mol:       the molecule of interest
      - tors_list: list of non-ring torsions
      - tors_list_rings: list of ring torsions
      - confId:    index of the conformation (default: first conformer)

      Return: list of torsion angles
  """
  torsions = []
  conf = mol.GetConformer(confId)
  for quartets, maxdev in tors_list:
    tors = []
    # loop over torsions and calculate angle
    for atoms in quartets:
      p1, p2, p3, p4 = _getTorsionAtomPositions(atoms, conf)
      tmpTors = (Geometry.ComputeSignedDihedralAngle(p1, p2, p3, p4) / math.pi) * 180.0
      if tmpTors < 0:
        tmpTors += 360.0  # angle between 0 and 360
      tors.append(tmpTors)
    torsions.append((tors, maxdev))
  # rings
  for quartets, maxdev in tors_list_rings:
    num = len(quartets)
    # loop over torsions and sum them up
    tors = 0
    for atoms in quartets:
      p1, p2, p3, p4 = _getTorsionAtomPositions(atoms, conf)
      tmpTors = abs((Geometry.ComputeSignedDihedralAngle(p1, p2, p3, p4) / math.pi) * 180.0)
      tors += tmpTors
    tors /= num
    torsions.append(([tors], maxdev))
  return torsions


def _findCentralBond(mol, distmat):
  """ Helper function to identify the atoms of the most central bond.

      Arguments:
      - mol:     the molecule of interest
      - distmat: distance matrix of the molecule

      Return: atom indices of the two most central atoms (in order)
  """
  from numpy import std

  # get the most central atom = atom with the least STD of shortest distances
  stds = []
  for i in range(mol.GetNumAtoms()):
    # only consider non-terminal atoms
    if len(_getHeavyAtomNeighbors(mol.GetAtomWithIdx(i))) < 2:
      continue
    tmp = [d for d in distmat[i]]
    tmp.pop(i)
    stds.append((std(tmp), i))
  stds.sort()
  aid1 = stds[0][1]
  # find the second most central bond that is bonded to aid1
  i = 1
  while 1:
    if mol.GetBondBetweenAtoms(aid1, stds[i][1]) is None:
      i += 1
    else:
      aid2 = stds[i][1]
      break
  return aid1, aid2  # most central atom comes first


def _calculateBeta(mol, distmat, aid1):
  """ Helper function to calculate the beta for torsion weights
      according to the formula in the paper.
      w(dmax/2) = 0.1

      Arguments:
      - mol:     the molecule of interest
      - distmat: distance matrix of the molecule
      - aid1:    atom index of the most central atom

      Return: value of beta (float)
  """
  # get all non-terminal bonds
  bonds = []
  for b in mol.GetBonds():
    nb1 = _getHeavyAtomNeighbors(b.GetBeginAtom())
    nb2 = _getHeavyAtomNeighbors(b.GetEndAtom())
    if len(nb2) > 1 and len(nb2) > 1:
      bonds.append(b)
  # get shortest distance
  dmax = 0
  for b in bonds:
    bid1 = b.GetBeginAtom().GetIdx()
    bid2 = b.GetEndAtom().GetIdx()
    d = max([distmat[aid1][bid1], distmat[aid1][bid2]])
    if (d > dmax):
      dmax = d
  dmax2 = dmax / 2.0
  beta = -math.log(0.1) / (dmax2 * dmax2)
  return beta


def CalculateTorsionWeights(mol, aid1=-1, aid2=-1, ignoreColinearBonds=True):
  """ Calculate the weights for the torsions in a molecule.
      By default, the highest weight is given to the bond 
      connecting the two most central atoms.
      If desired, two alternate atoms can be specified (must 
      be connected by a bond).

      Arguments:
      - mol:   the molecule of interest
      - aid1:  index of the first atom (default: most central)
      - aid2:  index of the second atom (default: second most central)
      - ignoreColinearBonds: if True (default), single bonds adjacent to
                             triple bonds are ignored
                             if False, alternative not-covalently bound
                             atoms are used to define the torsion

      Return: list of torsion weights (both non-ring and ring)
  """
  # get distance matrix
  distmat = Chem.GetDistanceMatrix(mol)
  if aid1 < 0 and aid2 < 0:
    aid1, aid2 = _findCentralBond(mol, distmat)
  else:
    b = mol.GetBondBetweenAtoms(aid1, aid2)
    if b is None:
      raise ValueError("Specified atoms must be connected by a bond.")
  # calculate beta according to the formula in the paper
  beta = _calculateBeta(mol, distmat, aid1)
  # get non-terminal, non-cyclic bonds
  bonds = _getBondsForTorsions(mol, ignoreColinearBonds)
  # get shortest paths and calculate weights
  weights = []
  for bid1, bid2, nb1, nb2 in bonds:
    if ((bid1, bid2) == (aid1, aid2)
        or (bid2, bid1) == (aid1, aid2)):  # if it's the most central bond itself
      d = 0
    else:
      # get shortest distance between the 4 atoms and add 1 to get bond distance
      d = min(distmat[aid1][bid1], distmat[aid1][bid2], distmat[aid2][bid1],
              distmat[aid2][bid2]) + 1
    w = math.exp(-beta * (d * d))
    weights.append(w)

  ## RINGS
  rings = mol.GetRingInfo()
  for r in rings.BondRings():
    # get shortest distances
    tmp = []
    num = len(r)
    for bidx in r:
      b = mol.GetBondWithIdx(bidx)
      bid1 = b.GetBeginAtomIdx()
      bid2 = b.GetEndAtomIdx()
      # get shortest distance between the 4 atoms and add 1 to get bond distance
      d = min(distmat[aid1][bid1], distmat[aid1][bid2], distmat[aid2][bid1],
              distmat[aid2][bid2]) + 1
      tmp.append(d)
    # calculate weights and append to list
    # Note: the description in the paper is not very clear, the following
    #       formula was found to give the same weights as shown in Fig. 1
    #       For a ring of size N: w = N/2 * exp(-beta*(sum(w of each bond in ring)/N)^2)
    w = sum(tmp) / float(num)
    w = math.exp(-beta * (w * w))
    weights.append(w * (num / 2.0))
  return weights


def CalculateTFD(torsions1, torsions2, weights=None):
  """ Calculate the torsion deviation fingerprint (TFD) given two lists of
      torsion angles.

      Arguments:
      - torsions1:  torsion angles of conformation 1
      - torsions2:  torsion angles of conformation 2
      - weights:    list of torsion weights (default: None)

      Return: TFD value (float)
  """
  if len(torsions1) != len(torsions2):
    raise ValueError("List of torsions angles must have the same size.")
  # calculate deviations and normalize (divide by max. possible deviation)
  deviations = []
  for tors1, tors2 in zip(torsions1, torsions2):
    mindiff = 180.0
    for t1 in tors1[0]:
      for t2 in tors2[0]:
        diff = abs(t1 - t2)
        if (360.0 - diff) < diff:  # we do not care about direction
          diff = 360.0 - diff
        #print t1, t2, diff
        if diff < mindiff:
          mindiff = diff
    deviations.append(mindiff / tors1[1])
  # do we use weights?
  if weights is not None:
    if len(weights) != len(torsions1):
      raise ValueError("List of torsions angles and weights must have the same size.")
    deviations = [d * w for d, w in zip(deviations, weights)]
    sum_weights = sum(weights)
  else:
    sum_weights = len(deviations)
  tfd = sum(deviations)
  if sum_weights != 0:  # avoid division by zero
    tfd /= sum_weights
  return tfd


def _getSameAtomOrder(mol1, mol2):
  """ Generate a new molecule with the atom order of mol1 and coordinates
      from mol2.
      
      Arguments:
      - mol1:     first instance of the molecule of interest
      - mol2:     second instance the molecule of interest

      Return: RDKit molecule
  """
  match = mol2.GetSubstructMatch(mol1)
  atomNums = tuple(range(mol1.GetNumAtoms()))
  if match != atomNums:  # atom orders are not the same!
    #print "Atoms of second molecule reordered."
    mol3 = Chem.Mol(mol1)
    mol3.RemoveAllConformers()
    for conf2 in mol2.GetConformers():
      confId = conf2.GetId()
      conf = rdchem.Conformer(mol1.GetNumAtoms())
      conf.SetId(confId)
      for i in range(mol1.GetNumAtoms()):
        conf.SetAtomPosition(i, mol2.GetConformer(confId).GetAtomPosition(match[i]))
      cid = mol3.AddConformer(conf)
    return mol3
  else:
    return Chem.Mol(mol2)


# some wrapper functions
def GetTFDBetweenConformers(mol, confIds1, confIds2, useWeights=True, maxDev='equal', symmRadius=2,
                            ignoreColinearBonds=True):
  """ Wrapper to calculate the TFD between two list of conformers 
      of a molecule

      Arguments:
      - mol:      the molecule of interest
      - confIds1:  first list of conformer indices
      - confIds2:  second list of conformer indices
      - useWeights: flag for using torsion weights in the TFD calculation
      - maxDev:   maximal deviation used for normalization
                  'equal': all torsions are normalized using 180.0 (default)
                  'spec':  each torsion is normalized using its specific
                           maximal deviation as given in the paper
      - symmRadius: radius used for calculating the atom invariants
                    (default: 2)
      - ignoreColinearBonds: if True (default), single bonds adjacent to
                             triple bonds are ignored
                             if False, alternative not-covalently bound
                             atoms are used to define the torsion

      Return: list of TFD values
  """
  tl, tlr = CalculateTorsionLists(mol, maxDev=maxDev, symmRadius=symmRadius,
                                  ignoreColinearBonds=ignoreColinearBonds)
  torsions1 = [CalculateTorsionAngles(mol, tl, tlr, confId=cid) for cid in confIds1]
  torsions2 = [CalculateTorsionAngles(mol, tl, tlr, confId=cid) for cid in confIds2]
  tfd = []
  if useWeights:
    weights = CalculateTorsionWeights(mol, ignoreColinearBonds=ignoreColinearBonds)
  else:
    weights = None
  for t1 in torsions1:
    for t2 in torsions2:
      tfd.append(CalculateTFD(t1, t2, weights=weights))
  return tfd


def GetTFDBetweenMolecules(mol1, mol2, confId1=-1, confId2=-1, useWeights=True, maxDev='equal',
                           symmRadius=2, ignoreColinearBonds=True):
  """ Wrapper to calculate the TFD between two molecules.
      Important: The two molecules must be isomorphic

      Arguments:
      - mol1:     first instance of the molecule of interest
      - mol2:     second instance the molecule of interest
      - confId1:  conformer index for mol1 (default: first conformer)
      - confId2:  conformer index for mol2 (default: first conformer)
      - useWeights: flag for using torsion weights in the TFD calculation
      - maxDev:   maximal deviation used for normalization
                  'equal': all torsions are normalized using 180.0 (default)
                  'spec':  each torsion is normalized using its specific
                           maximal deviation as given in the paper
      - symmRadius: radius used for calculating the atom invariants
                    (default: 2)
      - ignoreColinearBonds: if True (default), single bonds adjacent to
                             triple bonds are ignored
                             if False, alternative not-covalently bound
                             atoms are used to define the torsion

      Return: TFD value
  """
  if (Chem.MolToSmiles(mol1) != Chem.MolToSmiles(mol2)):
    raise ValueError("The two molecules must be instances of the same molecule!")
  mol2 = _getSameAtomOrder(mol1, mol2)
  tl, tlr = CalculateTorsionLists(mol1, maxDev=maxDev, symmRadius=symmRadius,
                                  ignoreColinearBonds=ignoreColinearBonds)
  # first molecule
  torsion1 = CalculateTorsionAngles(mol1, tl, tlr, confId=confId1)
  # second molecule
  torsion2 = CalculateTorsionAngles(mol2, tl, tlr, confId=confId2)
  if useWeights:
    weights = CalculateTorsionWeights(mol1, ignoreColinearBonds=ignoreColinearBonds)
  else:
    weights = None
  tfd = CalculateTFD(torsion1, torsion2, weights=weights)
  return tfd


def GetBestTFDBetweenMolecules(mol1, mol2, confId1=-1, useWeights=True, maxDev='equal',
                               symmRadius=2, ignoreColinearBonds=True):
  """ Wrapper to calculate the best TFD between a single conformer of mol1 and all the conformers of mol2
      Important: The two molecules must be isomorphic

      Arguments:
      - mol1:     first instance of the molecule of interest
      - mol2:     second instance the molecule of interest
      - confId1:  conformer index for mol1 (default: first conformer)
      - useWeights: flag for using torsion weights in the TFD calculation
      - maxDev:   maximal deviation used for normalization
                  'equal': all torsions are normalized using 180.0 (default)
                  'spec':  each torsion is normalized using its specific
                           maximal deviation as given in the paper
      - symmRadius: radius used for calculating the atom invariants
                    (default: 2)
      - ignoreColinearBonds: if True (default), single bonds adjacent to
                             triple bonds are ignored
                             if False, alternative not-covalently bound
                             atoms are used to define the torsion

      Return: TFD value
  """
  if (Chem.MolToSmiles(mol1) != Chem.MolToSmiles(mol2)):
    raise ValueError("The two molecules must be instances of the same molecule!")
  mol2 = _getSameAtomOrder(mol1, mol2)
  tl, tlr = CalculateTorsionLists(mol1, maxDev=maxDev, symmRadius=symmRadius,
                                  ignoreColinearBonds=ignoreColinearBonds)
  # first molecule
  torsion1 = CalculateTorsionAngles(mol1, tl, tlr, confId=confId1)
  if useWeights:
    weights = CalculateTorsionWeights(mol1, ignoreColinearBonds=ignoreColinearBonds)
  else:
    weights = None
  best = 1e8
  for conf in mol2.GetConformers():
    # second molecule
    torsion2 = CalculateTorsionAngles(mol2, tl, tlr, confId=conf.GetId())
    tfd = CalculateTFD(torsion1, torsion2, weights=weights)
    best = min(best, tfd)
  return best


def GetTFDMatrix(mol, useWeights=True, maxDev='equal', symmRadius=2, ignoreColinearBonds=True):
  """ Wrapper to calculate the matrix of TFD values for the
      conformers of a molecule.

      Arguments:
      - mol:      the molecule of interest
      - useWeights: flag for using torsion weights in the TFD calculation
      - maxDev:   maximal deviation used for normalization
                  'equal': all torsions are normalized using 180.0 (default)
                  'spec':  each torsion is normalized using its specific
                           maximal deviation as given in the paper
      - symmRadius: radius used for calculating the atom invariants
                    (default: 2)
      - ignoreColinearBonds: if True (default), single bonds adjacent to
                             triple bonds are ignored
                             if False, alternative not-covalently bound
                             atoms are used to define the torsion

      Return: matrix of TFD values
      Note that the returned matrix is symmetrical, i.e. it is the
      lower half of the matrix, e.g. for 5 conformers:
      matrix = [ a,
                 b, c,
                 d, e, f,
                 g, h, i, j]
  """
  tl, tlr = CalculateTorsionLists(mol, maxDev=maxDev, symmRadius=symmRadius,
                                  ignoreColinearBonds=ignoreColinearBonds)
  numconf = mol.GetNumConformers()
  torsions = [
    CalculateTorsionAngles(mol, tl, tlr, confId=conf.GetId()) for conf in mol.GetConformers()
  ]
  tfdmat = []
  if useWeights:
    weights = CalculateTorsionWeights(mol, ignoreColinearBonds=ignoreColinearBonds)
  else:
    weights = None
  for i in range(0, numconf):
    for j in range(0, i):
      tfdmat.append(CalculateTFD(torsions[i], torsions[j], weights=weights))

  return tfdmat
