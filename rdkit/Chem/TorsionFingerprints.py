# $Id$
#
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" Torsion Fingerprints

"""
from rdkit import rdBase
from rdkit import RDConfig
from rdkit import Geometry
from rdkit import Chem
from rdkit.Chem import AllChem
import math, os, copy
from rdkit.RDLogger import logger
logger = logger()
import warnings

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
    fp = AllChem.GetMorganFingerprint(mol, radius, fromAtoms=[i], bitInfo=info)
    for k in info.keys():
      if info[k][0][1] == radius:
        inv.append(k)
  return inv

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
                  'equal': all torsions are normalized using 180.0
                  'spec':  each torsion is normalized using its specific
                           maximal deviation as given in the paper
      - symmRadius: radius used for calculating the atom invariants

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
    nb1 = [n for n in b.GetBeginAtom().GetNeighbors() if (n.GetSymbol()!='H' and n.GetIdx()!=a2)]
    nb2 = [n for n in b.GetEndAtom().GetNeighbors() if (n.GetSymbol()!='H' and n.GetIdx()!=a1)]
    if nb1 and nb2: # no terminal bonds
      bonds.append((b, a1, a2, nb1, nb2))
  # get atom invariants
  if symmRadius > 0:
    inv = _getAtomInvariantsWithRadius(mol, symmRadius)
  else:
    inv = AllChem.GetConnectivityInvariants(mol)
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
      if len(nb1) == 2 and len(nb2) == 2: # case 9
        tors_list.append(([[(n1.GetIdx(), a1, a2, n2.GetIdx()) for n1 in d1] for n2 in d2], 90.0))
      elif len(nb1) == 3 and len(nb2) == 3: # case 21
        tors_list.append(([[(n1.GetIdx(), a1, a2, n2.GetIdx()) for n1 in d1] for n2 in d2], 60.0))
      else: # case 15
        tors_list.append(([[(n1.GetIdx(), a1, a2, n2.GetIdx()) for n1 in d1] for n2 in d2], 30.0))
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
  """ Calculate the torsion angels for a list of non-ring and 
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
    else:
      # loop over torsions and take minimum
      tors = 360.0
      for t2 in t:
        p1, p2, p3, p4 = _getTorsionAtomPositions(t2, conf)
        tmp = (Geometry.ComputeSignedDihedralAngle(p1, p2, p3, p4)/math.pi)*180.0
        if abs(tmp) < tors: tors = abs(tmp)
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
    tmp = [d for d in distmat[i]]
    tmp.pop(i)
    stds.append((std(tmp), i))
  stds.sort()
  aid1 = stds[0][1]
  # find the second most central bond that is bonded to aid1
  nbs = [nb.GetIdx() for nb in [n for n in mol.GetAtomWithIdx(aid1).GetNeighbors() if n.GetSymbol()!='H'] if len([n for n in nb.GetNeighbors() if n.GetSymbol()!='H']) > 0]
  nbs = [s for s in stds if s[1] in nbs]
  aid2 = nbs[0][1]
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
    nb1 = [n for n in b.GetBeginAtom().GetNeighbors() if n.GetSymbol()!='H']
    nb2 = [n for n in b.GetEndAtom().GetNeighbors() if n.GetSymbol()!='H']
    if nb2 > 1 and nb2 > 1:
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
  distmat = AllChem.GetDistanceMatrix(mol)
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
    nb1 = [n for n in b.GetBeginAtom().GetNeighbors() if n.GetSymbol()!='H']
    nb2 = [n for n in b.GetEndAtom().GetNeighbors() if n.GetSymbol()!='H']
    if nb1 > 1 and nb2 > 1:
      bonds.append(b)
  # get shortest paths
  shortest_distances = []
  for b in bonds:
    bid1 = b.GetBeginAtom().GetIdx()
    bid2 = b.GetEndAtom().GetIdx()
    if (set([bid1, bid2]) == set([aid1, aid2])):
      d = 0
    else:
      # get shortest distance between the 4 atoms and add 1 to get bond distance
      d = {distmat[aid1][bid1], distmat[aid1][bid2], distmat[aid2][bid1], distmat[aid2][bid2]}.pop()+1
    shortest_distances.append(d)
  # calculate weights
  weights = []
  for d in shortest_distances:
    w = math.exp(-beta*(d*d))
    weights.append(w)

  ## RINGS
  shortest_distances = []
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
      d = {distmat[aid1][bid1], distmat[aid1][bid2], distmat[aid2][bid1], distmat[aid2][bid2]}.pop()+1
      tmp.append(d)
    shortest_distances.append(tmp)
  # calculate weights and append to list
  for dists in shortest_distances:
    w = [math.exp(-beta*(d*d)) for d in dists]
    w = sum(w)/2.0
    weights.append(w)
  return weights

def CalculateTDF(torsions1, torsions2, weights=None):
  """ Calculate the torsion deviation fingerprint (TDF) given two lists of
      torsion angles.

      Arguments;
      - torsions1:  torsion angles of conformation 1
      - torsions2:  torsion angles of conformation 2
      - weights:    list of torsion weights (default: None)

      Return: TDF value (float)
  """
  if len(torsions1) != len(torsions2):
    raise ValueError("List of torsions angles must have the same size.")
  # calculate deviations and normalize (divide by max. possible deviation)
  deviations = []
  for t1, t2 in zip(torsions1, torsions2):
    deviations.append(abs(t1[0]-t2[0])/t1[1])
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
