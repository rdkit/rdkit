#  Copyright (C) 2014 Sereina Riniker
#
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" Torsion Fingerprints (Deviation) (TFD)
    According to a paper from Schulz-Gasch et al., JCIM, 52, 1499-1512 (2012).

"""
from rdkit import rdBase
from rdkit import RDConfig
from rdkit import Geometry
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import math, os

def _doMatch(inv, atoms):
  """ Helper function to check if all atoms in the list are the same
      
      Arguments:
      - inv:    atom invariants (used to define equivalence of atoms)
      - atoms:  list of atoms to check

      Return: boolean
  """
  match = True
  for i in range(len(atoms)-1):
    for j in range(i+1, len(atoms)):
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
  for i in range(len(atoms)-1):
    for j in range(i+1, len(atoms)):
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
  for i in range(mol.GetNumAtoms()):
    info = {}
    fp = rdMolDescriptors.GetMorganFingerprint(mol, radius, fromAtoms=[i], bitInfo=info)
    for k in info.keys():
      if info[k][0][1] == radius:
        inv.append(k)
  return inv

def _getHeavyAtomNeighbors(atom1, aid2=-1):
  """ Helper function to calculate the number of heavy atom neighbors.

      Arguments:
      - atom1:    the atom of interest
      - aid2:     atom index that should be excluded from neighbors (default: none)

      Return: a list of heavy atom neighbors of the given atom
  """
  if aid2 < 0:
    return [n for n in atom1.GetNeighbors() if n.GetSymbol()!='H']
  else:
    return [n for n in atom1.GetNeighbors() if (n.GetSymbol()!='H' and n.GetIdx()!=aid2)]

def _getIndexforTorsion(neighbors, inv):
  """ Helper function to calculate the index of the reference atom for 
      a given atom

      Arguments:
      - neighbors:  list of the neighbors of the atom
      - inv:        atom invariants

      Return: list of atom indices as reference for torsion
  """
  if len(neighbors) == 1: # atom has only one neighbor
    return [neighbors[0]]
  elif _doMatch(inv, neighbors): # atom has all symmetric neighbors
    return neighbors
  elif _doNotMatch(inv, neighbors): # atom has all different neighbors
    # simply use the first neighbor
    return [neighbors[0]]
  at = _doMatchExcept1(inv, neighbors) # two neighbors the same, one different
  if at is None:
    raise ValueError("Atom neighbors are either all the same or all different")
  return [at] 

def CalculateTorsionLists(mol, maxDev='equal', symmRadius=2):
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

      Return: two lists of torsions: non-ring and ring torsions
  """
  if maxDev not in ['equal', 'spec']:
    raise ValueError("maxDev must be either equal or spec")
  # get non-terminal, non-cyclic bonds
  bonds = []
  for b in mol.GetBonds():
    if b.IsInRing(): continue
    a1 = b.GetBeginAtomIdx()
    a2 = b.GetEndAtomIdx()
    nb1 = _getHeavyAtomNeighbors(b.GetBeginAtom(), a2)
    nb2 = _getHeavyAtomNeighbors(b.GetEndAtom(), a1)
    if nb1 and nb2: # no terminal bonds
      bonds.append((b, a1, a2, nb1, nb2))
  # get atom invariants
  if symmRadius > 0:
    inv = _getAtomInvariantsWithRadius(mol, symmRadius)
  else:
    inv = rdMolDescriptors.GetConnectivityInvariants(mol)
  # get the torsions
  tors_list = [] # to store the atom indices of the torsions
  for b, a1, a2, nb1, nb2 in bonds:
    d1 = _getIndexforTorsion(nb1, inv)
    d2 = _getIndexforTorsion(nb2, inv)
    if len(d1) == 1 and len(d2) == 1: # case 1, 2, 4, 5, 7, 10, 16, 12, 17, 19
      tors_list.append(([(d1[0].GetIdx(), a1, a2, d2[0].GetIdx())], 180.0))
    elif len(d1) == 1: # case 3, 6, 8, 13, 20
      if len(nb2) == 2: # two neighbors
        tors_list.append(([(d1[0].GetIdx(), a1, a2, nb.GetIdx()) for nb in d2], 90.0))
      else: # three neighbors
        tors_list.append(([(d1[0].GetIdx(), a1, a2, nb.GetIdx()) for nb in d2], 60.0))
    elif len(d2) == 1: # case 3, 6, 8, 13, 20
      if len(nb1) == 2:
        tors_list.append(([(nb.GetIdx(), a1, a2, d2[0].GetIdx()) for nb in d1], 90.0))
      else: # three neighbors
        tors_list.append(([(nb.GetIdx(), a1, a2, d2[0].GetIdx()) for nb in d1], 60.0))
    else: # both symmetric
      tmp = []
      for n1 in d1:
        for n2 in d2:
          tmp.append((n1.GetIdx(), a1, a2, n2.GetIdx()))
      if len(nb1) == 2 and len(nb2) == 2: # case 9
        tors_list.append((tmp, 90.0))
      elif len(nb1) == 3 and len(nb2) == 3: # case 21
        tors_list.append((tmp, 60.0))
      else: # case 15
        tors_list.append((tmp, 30.0))
  # maximal possible deviation for non-cyclic bonds
  if maxDev == 'equal':
    tors_list = [(t,180.0) for t,d in tors_list] 
  # rings
  rings = Chem.GetSymmSSSR(mol)
  tors_list_rings = []
  for r in rings:
    # get the torsions
    tmp = []
    num = len(r)
    maxdev = 180.0 * math.exp(-0.025*(num-14)*(num-14))
    for i in range(len(r)):
      tmp.append((r[i], r[(i+1)%num], r[(i+2)%num], r[(i+3)%num]))
    tors_list_rings.append((tmp,maxdev))
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
  for t,maxdev in tors_list:
    if len(t) == 1:
      t = t[0]
      p1, p2, p3, p4 = _getTorsionAtomPositions(t, conf)
      tors = (Geometry.ComputeSignedDihedralAngle(p1, p2, p3, p4)/math.pi)*180.0
      if tors < 0: tors += 360.0 # angle between 0 and 360
    else:
      # loop over torsions and take minimum
      tors = 360.0
      for t2 in t:
        p1, p2, p3, p4 = _getTorsionAtomPositions(t2, conf)
        tmp = (Geometry.ComputeSignedDihedralAngle(p1, p2, p3, p4)/math.pi)*180.0
        if tmp < 0: tmp += 360.0 # angle between 0 and 360
        if tmp < tors: tors = tmp
    torsions.append((tors, maxdev))
  # rings
  for t,maxdev in tors_list_rings:
    num = len(t)
    # loop over torsions and sum them up
    tors = 0
    for t2 in t:
      p1, p2, p3, p4 = _getTorsionAtomPositions(t2, conf)
      tmp = abs((Geometry.ComputeSignedDihedralAngle(p1, p2, p3, p4)/math.pi)*180.0)
      tors += tmp
    tors /= num
    torsions.append((tors, maxdev))
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
    if len(_getHeavyAtomNeighbors(mol.GetAtomWithIdx(i))) < 2: continue
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
  return aid1, aid2 # most central atom comes first

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
    if (d > dmax): dmax = d
  dmax2 = dmax/2.0
  beta = -math.log(0.1)/(dmax2*dmax2)
  return beta

def CalculateTorsionWeights(mol, aid1=-1, aid2=-1):
  """ Calculate the weights for the torsions in a molecule.
      By default, the highest weight is given to the bond 
      connecting the two most central atoms.
      If desired, two alternate atoms can be specified (must 
      be connected by a bond).

      Arguments:
      - mol:   the molecule of interest
      - aid1:  index of the first atom (default: most central)
      - aid2:  index of the second atom (default: second most cenral)

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
  bonds = []
  for b in mol.GetBonds():
    if b.IsInRing(): continue
    nb1 = _getHeavyAtomNeighbors(b.GetBeginAtom())
    nb2 = _getHeavyAtomNeighbors(b.GetEndAtom())
    if len(nb1) > 1 and len(nb2) > 1:
      bonds.append(b)
  # get shortest paths and calculate weights
  weights = []
  for b in bonds:
    bid1 = b.GetBeginAtom().GetIdx()
    bid2 = b.GetEndAtom().GetIdx()
    if ((bid1, bid2) == (aid1, aid2)
      or (bid2, bid1) == (aid1, aid2)): # if it's the most central bond itself
      d = 0
    else:
      # get shortest distance between the 4 atoms and add 1 to get bond distance
      d = min(distmat[aid1][bid1], distmat[aid1][bid2], distmat[aid2][bid1], distmat[aid2][bid2])+1
    w = math.exp(-beta*(d*d))
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
      d = min(distmat[aid1][bid1], distmat[aid1][bid2], distmat[aid2][bid1], distmat[aid2][bid2])+1
      tmp.append(d)
    # calculate weights and append to list
    # Note: the description in the paper is not very clear, the following
    #       formula was found to give the same weights as shown in Fig. 1
    #       For a ring of size N: w = N/2 * exp(-beta*(sum(w of each bond in ring)/N)^2)
    w = sum(tmp)/float(num)
    w = math.exp(-beta*(w*w))
    weights.append(w*(num/2.0))
  return weights

def CalculateTFD(torsions1, torsions2, weights=None):
  """ Calculate the torsion deviation fingerprint (TFD) given two lists of
      torsion angles.

      Arguments;
      - torsions1:  torsion angles of conformation 1
      - torsions2:  torsion angles of conformation 2
      - weights:    list of torsion weights (default: None)

      Return: TFD value (float)
  """
  if len(torsions1) != len(torsions2):
    raise ValueError("List of torsions angles must have the same size.")
  # calculate deviations and normalize (divide by max. possible deviation)
  deviations = []
  for t1, t2 in zip(torsions1, torsions2):
    diff = abs(t1[0]-t2[0])
    if (360.0-diff) < diff: # we do not care about direction
      diff = 360.0 - diff
    deviations.append(diff/t1[1])
  # do we use weights?
  if weights is not None:
    if len(weights) != len(torsions1):
      raise ValueError("List of torsions angles and weights must have the same size.")
    deviations = [d*w for d,w in zip(deviations, weights)]
    sum_weights = sum(weights)
  else:
    sum_weights = len(deviations)
  tfd = sum(deviations)
  if sum_weights != 0: # avoid division by zero
    tfd /= sum_weights
  return tfd

# some wrapper functions
def GetTFDBetweenConformers(mol, confIds1, confIds2, useWeights=True, maxDev='equal', symmRadius=2):
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

      Return: list of TFD values
  """
  tl, tlr = CalculateTorsionLists(mol, maxDev=maxDev, symmRadius=symmRadius)
  torsions1 = [CalculateTorsionAngles(mol, tl, tlr, confId=cid) for cid in confIds1]
  torsions2 = [CalculateTorsionAngles(mol, tl, tlr, confId=cid) for cid in confIds2]
  tfd = []
  if useWeights:
    weights = CalculateTorsionWeights(mol)
    for t1 in torsions1:
      for t2 in torsions2:
        tfd.append(CalculateTFD(t1, t2, weights=weights))
  else:
    for t1 in torsions1:
      for t2 in torsions2:
        tfd.append(CalculateTFD(t1, t2))
  return tfd

def GetTFDBetweenMolecules(mol1, mol2, confIds1=-1, confIds2=-1, useWeights=True, maxDev='equal', symmRadius=2):
  """ Wrapper to calculate the TFD between two list of conformers 
      of two molecules.
      Important: The two molecules must be instances of the same molecule

      Arguments:
      - mol1:     first instance of the molecule of interest
      - mol2:     second instance the molecule of interest
      - confIds1:  list of conformer indices from mol1 (default: first conformer)
      - confIds2:  list of conformer indices from mol2 (default: first conformer)
      - useWeights: flag for using torsion weights in the TFD calculation
      - maxDev:   maximal deviation used for normalization
                  'equal': all torsions are normalized using 180.0 (default)
                  'spec':  each torsion is normalized using its specific
                           maximal deviation as given in the paper
      - symmRadius: radius used for calculating the atom invariants
                    (default: 2)

      Return: list of TFD values
  """
  if (Chem.MolToSmiles(mol1) != Chem.MolToSmiles(mol2)):
    raise ValueError("The two molecules must be instances of the same molecule!")
  tl, tlr = CalculateTorsionLists(mol1, maxDev=maxDev, symmRadius=symmRadius)
  # first molecule
  if confIds1 < 0:
    torsions1 = [CalculateTorsionAngles(mol1, tl, tlr)]
  else:
    torsions1 = [CalculateTorsionAngles(mol1, tl, tlr, confId=cid) for cid in confIds1]
  # second molecule
  if confIds2 < 0:
    torsions2 = [CalculateTorsionAngles(mol2, tl, tlr)]
  else:
    torsions2 = [CalculateTorsionAngles(mol2, tl, tlr, confId=cid) for cid in confIds2]
  tfd = []
  if useWeights:
    weights = CalculateTorsionWeights(mol1)
    for t1 in torsions1:
      for t2 in torsions2:
        tfd.append(CalculateTFD(t1, t2, weights=weights))
  else:
    for t1 in torsions1:
      for t2 in torsions2:
        tfd.append(CalculateTFD(t1, t2))
  return tfd

def GetTFDMatrix(mol, useWeights=True, maxDev='equal', symmRadius=2):
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

      Return: matrix of TFD values
      Note that the returned matrix is symmetrical, i.e. it is the
      lower half of the matrix, e.g. for 5 conformers:
      matrix = [ a,
                 b, c,
                 d, e, f,
                 g, h, i, j]
  """
  tl, tlr = CalculateTorsionLists(mol, maxDev=maxDev, symmRadius=symmRadius)
  numconf = mol.GetNumConformers()
  torsions = [CalculateTorsionAngles(mol, tl, tlr, confId=conf.GetId()) for conf in mol.GetConformers()]
  tfdmat = []
  if useWeights:
    weights = CalculateTorsionWeights(mol)
    for i in range(0, numconf):
      for j in range(0, i):
        tfdmat.append(CalculateTFD(torsions[i], torsions[j], weights=weights))
  else:
    for i in range(0, numconf):
      for j in range(0, i):
        tfdmat.append(CalculateTFD(torsions[i], torsions[j]))
  return tfdmat

