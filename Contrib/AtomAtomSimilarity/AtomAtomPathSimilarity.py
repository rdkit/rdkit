# This is a reference implementation of the Atom-Atom-Path (AAP) similarity metric
#  from Gobbi et al, J. Cheminf. (2015) 7:11. https://doi.org/10.1186/s13321-015-0056-8
#
# Original author: Richard Hall
#
import time
import unittest

import numpy
from scipy.optimize import linear_sum_assignment

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, rdmolops
from rdkit.Chem.Fingerprints import FingerprintMols

_BK_ = {
  Chem.rdchem.BondType.SINGLE: 1,
  Chem.rdchem.BondType.DOUBLE: 2,
  Chem.rdchem.BondType.TRIPLE: 3,
  Chem.rdchem.BondType.AROMATIC: 4
}
_BONDSYMBOL_ = {1: '-', 2: '=', 3: '#', 4: ':'}

#_nAT_ = 217 # 108*2+1
_nAT_ = 223  # Gobbi code actually uses the first prime higher than 217, not 217 itself
_nBT_ = 5

#def FindAllPathsOfLengthN_Gobbi(mol, length, rootedAtAtom=-1, uniquepaths=True):
#	return FindAllPathsOfLengthMToN(mol, length, length, rootedAtAtom=rootedAtAtom, uniquepaths=uniquepaths)


def FindAllPathsOfLengthMToN_Gobbi(mol, minlength, maxlength, rootedAtAtom=-1, uniquepaths=True):
  '''this function returns the same set of bond paths as the Gobbi paper.  These differ a little from the rdkit FindAllPathsOfLengthMToN function'''
  paths = []
  for atom in mol.GetAtoms():
    if rootedAtAtom == -1 or atom.GetIdx() == rootedAtAtom:
      path = []
      visited = set([atom.GetIdx()])
      #			visited = set()
      _FindAllPathsOfLengthMToN_Gobbi(atom, path, minlength, maxlength, visited, paths)

  if uniquepaths:
    uniquepathlist = []
    seen = set()
    for path in paths:
      if path not in seen:
        reversepath = tuple([i for i in path[::-1]])
        if reversepath not in seen:
          uniquepathlist.append(path)
          seen.add(path)
    return uniquepathlist
  else:
    return paths


def _FindAllPathsOfLengthMToN_Gobbi(atom, path, minlength, maxlength, visited, paths):
  for bond in atom.GetBonds():
    if bond.GetIdx() not in path:
      bidx = bond.GetIdx()
      path.append(bidx)
      if len(path) >= minlength and len(path) <= maxlength:
        paths.append(tuple(path))
      if len(path) < maxlength:
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetIdx() == atom.GetIdx():
          nextatom = a2
        else:
          nextatom = a1
        nextatomidx = nextatom.GetIdx()
        if nextatomidx not in visited:
          visited.add(nextatomidx)
          _FindAllPathsOfLengthMToN_Gobbi(nextatom, path, minlength, maxlength, visited, paths)
          visited.remove(nextatomidx)
      path.pop()


def getpathintegers(m1, uptolength=7):
  '''returns a list of integers describing the paths for molecule m1.  This uses numpy 16 bit unsigned integers to reproduce the data in the Gobbi paper.  The returned list is sorted'''
  bondtypelookup = {}
  for b in m1.GetBonds():
    bondtypelookup[b.GetIdx()] = _BK_[b.GetBondType()], b.GetBeginAtom(), b.GetEndAtom()
  pathintegers = {}
  for a in m1.GetAtoms():
    idx = a.GetIdx()
    pathintegers[idx] = []
    #		for pathlength in range(1, uptolength+1):
    #			for path in rdmolops.FindAllPathsOfLengthN(m1, pathlength, rootedAtAtom=idx):
    for ipath, path in enumerate(
        FindAllPathsOfLengthMToN_Gobbi(m1, 1, uptolength, rootedAtAtom=idx, uniquepaths=False)):
      strpath = []
      currentidx = idx
      res = []
      for ip, p in enumerate(path):
        bk, a1, a2 = bondtypelookup[p]
        strpath.append(_BONDSYMBOL_[bk])
        if a1.GetIdx() == currentidx:
          a = a2
        else:
          a = a1
        ak = a.GetAtomicNum()
        if a.GetIsAromatic():
          ak += 108
#trying to get the same behaviour as the Gobbi test code - it looks like a circular path includes the bond, but not the closure atom - this fix works
        if a.GetIdx() == idx:
          ak = None
        if ak is not None:
          astr = a.GetSymbol()
          if a.GetIsAromatic():
            strpath.append(astr.lower())
          else:
            strpath.append(astr)
        res.append((bk, ak))
        currentidx = a.GetIdx()
      pathuniqueint = numpy.ushort(0)  # work with 16 bit unsigned integers and ignore overflow...
      for ires, (bi, ai) in enumerate(res):
        #use 16 bit unsigned integer arithmetic to reproduce the Gobbi ints
        #					pathuniqueint = ((pathuniqueint+bi)*_nAT_+ai)*_nBT_
        val1 = pathuniqueint + numpy.ushort(bi)
        val2 = val1 * numpy.ushort(_nAT_)
        #trying to get the same behaviour as the Gobbi test code - it looks like a circular path includes the bond, but not the closure atom - this fix works
        if ai is not None:
          val3 = val2 + numpy.ushort(ai)
          val4 = val3 * numpy.ushort(_nBT_)
        else:
          val4 = val2
        pathuniqueint = val4
      pathintegers[idx].append(pathuniqueint)


#sorted lists allow for a quicker comparison algorithm
  for p in pathintegers.values():
    p.sort()
  return pathintegers


def getcommon(l1, ll1, l2, ll2):
  '''returns the number of items sorted lists l1 and l2 have in common.  ll1 and ll2 are the list lengths'''
  ncommon = 0
  ix1 = 0
  ix2 = 0
  while (ix1 < ll1) and (ix2 < ll2):
    a1 = l1[ix1]
    a2 = l2[ix2]
    #a1 is < or > more often that ==
    if a1 < a2:
      ix1 += 1
    elif a1 > a2:
      ix2 += 1
    else:  # a1 == a2:
      ncommon += 1
      ix1 += 1
      ix2 += 1
  return ncommon


def getsimaibj(aipaths, bjpaths, naipaths, nbjpaths):
  '''returns the similarity of two sorted path lists.  Equation 2'''
  nc = getcommon(aipaths, naipaths, bjpaths, nbjpaths)
  sim = float(nc + 1) / (max(naipaths, nbjpaths) * 2 - nc + 1)
  return sim


def getmappings(simmatrixarray):
  '''return a mapping of the atoms in the similarity matix using the heuristic algorithm described in the paper'''

  costarray = numpy.ones(simmatrixarray.shape) - simmatrixarray

  it = numpy.nditer(costarray, flags=['multi_index'], op_flags=['writeonly'])
  dsu = []
  for a in it:
    dsu.append((a, it.multi_index[0], it.multi_index[1]))
  dsu.sort()

  seena = set()
  seenb = set()
  mappings = []
  for sim, a, b in dsu:
    if a not in seena and b not in seenb:
      seena.add(a)
      seenb.add(b)
      mappings.append((a, b))

  return mappings[:min(simmatrixarray.shape)]


def gethungarianmappings(simmatrixarray):
  '''return a mapping of the atoms in the similarity matrix - the Hungarian algorithm is used because it is invariant to atom ordering.  Requires scipy'''
  costarray = numpy.ones(simmatrixarray.shape) - simmatrixarray
  row_ind, col_ind = linear_sum_assignment(costarray)
  res = zip(row_ind, col_ind)
  return res


def getsimab(mappings, simmatrixdict):
  '''return the similarity for a set of mapping.  See Eqn 3'''
  naa, nab = simmatrixdict.shape

  score = 0.0
  for a, b in mappings:
    score += simmatrixdict[a][b]
  simab = score / (max(naa, nab) * 2 - score)
  return simab


def getsimmatrix(m1, m1pathintegers, m2, m2pathintegers):
  '''generate a matrix of atom atom similarities.  See Figure 4'''

  aidata = [((ai.GetAtomicNum(), ai.GetIsAromatic()), ai.GetIdx()) for ai in m1.GetAtoms()]
  bjdata = [((bj.GetAtomicNum(), bj.GetIsAromatic()), bj.GetIdx()) for bj in m2.GetAtoms()]

  simmatrixarray = numpy.zeros((len(aidata), len(bjdata)))

  for ai, (aitype, aiidx) in enumerate(aidata):
    aipaths = m1pathintegers[aiidx]
    naipaths = len(aipaths)
    for bj, (bjtype, bjidx) in enumerate(bjdata):
      if aitype == bjtype:
        bjpaths = m2pathintegers[bjidx]
        nbjpaths = len(bjpaths)
        simmatrixarray[ai][bj] = getsimaibj(aipaths, bjpaths, naipaths, nbjpaths)
  return simmatrixarray


def AtomAtomPathSimilarity(m1, m2, m1pathintegers=None, m2pathintegers=None):
  '''compute the Atom Atom Path Similarity for a pair of RDKit molecules.  See Gobbi et al, J. ChemInf (2015) 7:11
      the most expensive part of the calculation is computing the path integers - we can precompute these and pass them in as an argument'''
  if m1pathintegers is None:
    m1pathintegers = getpathintegers(m1)
  if m2pathintegers is None:
    m2pathintegers = getpathintegers(m2)

  simmatrix = getsimmatrix(m1, m1pathintegers, m2, m2pathintegers)

  #	mappings = getmappings(simmatrix)
  mappings = gethungarianmappings(simmatrix)

  simab = getsimab(mappings, simmatrix)

  return simab


def test0():
  '''reproduce the worked similarity in the Gobbi paper'''
  m1 = Chem.MolFromSmiles('o1nccc1C')
  m2 = Chem.MolFromSmiles('[nH]1nccc1')
  return AtomAtomPathSimilarity(m1, m2)


def test1():
  '''generate a set of path integers for 3 molecules from the Gobbi source IAAPathGeneratorCharTest.java'''
  res = []
  smiles = ["C", "C(=O)F", "C1ON1"]
  for s in smiles:
    m = Chem.MolFromSmiles(s)
    mpathintegers = getpathintegers(m)
    res.append(mpathintegers)
  return res


def test2():
  '''generate a matrix molecules from the Gobbi source AAPathComparator2Test.java'''
  smileslist = [
    "*", "C", "N", "CCO", "CC(=O)N", "c1ccccc1", "c1ncncc1", "c1[nH]ccc1", "c1ncncc1CC(=O)N",
    "c1ccccc1c1ncncc1"
  ]
  sims = []
  for s1 in smileslist:
    for s2 in smileslist:
      m1 = Chem.MolFromSmiles(s1)
      m2 = Chem.MolFromSmiles(s2)
      sims.append('%.4f' % AtomAtomPathSimilarity(m1, m2))
  return sims


def test3():
  '''generate a set of similarities for the example compounds in Figure 1.  These are compared to the values in Additional File 1'''
  m1a = Chem.MolFromSmiles('Clc1ccc(CN2CCC(CC2)c3cc([nH]n3)c4ccc(Cl)cc4)cc1')
  m1b = Chem.MolFromSmiles('Clc1ccc(CN2CCN(CC2)CC(=O)N(C)c3ccccc3)cc1')
  m2a = Chem.MolFromSmiles('Cc1cccn2cc(nc12)c3ccc(NC(=O)CN4CCCC4)cc3')
  m2b = Chem.MolFromSmiles('Cc1c(cc2ccccn12)c3ccc(OCCCN4CCCCC4)cc3')

  res = []
  for m1 in (m1a, m1b, m2a, m2b):
    for m2 in (m1a, m1b, m2a, m2b):
      sim = AtomAtomPathSimilarity(m1, m2)
      res.append('%.2f' % sim)
  return res


def timeit():
  #these are the first 40 smiles in the Gobbi 13321_2015_56_MOESM2_ESM file'''
  molstr = '''C[C@@H](O)[C@@H]1OCC[C@@H](C)[C@H](O)C(=O)OC[C@]23CCC(C)=C[C@H]2O[C@@H]4C[C@@H](OC(=O)C=CC=C1)[C@@]3(C)[C@@]45CO5
CC1=C[C@H]2O[C@@H]3C[C@H]4OC(=O)C=CC=CC(=O)OCCC(C)=CC(=O)OC[C@@]2(CC1)[C@]4(C)[C@]35CO5
CC1=CC(=O)OC[C@]23C[C@H](O)C(C)=C[C@H]2O[C@@H]4C[C@@H](OC(=O)C=CC=CC(=O)OCC1)[C@@]3(C)[C@]45CO5
CC1(C)N=C(N)N=C(N)N1C2=CC=C(Br)C=C2
CC1(C)N=C(N)N=C(N)N1C2=CC=CC=C2
CC1(C)N=C(N)N=C(N)N1C2=CC=C(I)C=C2
CC1=CC=C(C=C1)N2C(N)=NC(N)=NC2(C)C
CC1(C)N=C(N)N=C(N)N1C2=CC=C(Cl)C=C2
CC1(C)N=C(N)N=C(N)N1C2=CC=C(F)C=C2
CC1=CC=CC(=C1)N2C(N)=NC(N)=NC2(C)C
COC1=CC=C(C=C1)N2C(N)=NC(N)=NC2(C)C
CC1=CC=C(N2C(N)=NC(N)=NC2(C)C)C(C)=C1
CCOC1=CC=C(C=C1)N2C(N)=NC(N)=NC2(C)C
COC1=CC=CC(=C1)N2C(N)=NC(N)=NC2(C)C
CC1=CC=CC(NC2=NC(N)=NC(C)(C)N2)=C1
CNC1=C(N(CC2=CC=C(Cl)C(Cl)=C2)C(C)=O)C(=O)C3=CC=CC=C3C1=O
CC(=O)N(CC1=CC=C(F)C=C1)C2=C(NCC3=CC=CC=C3)C(=O)C4=CC=CC=C4C2=O
CC(=O)N(CC1=CC=C(F)C=C1)C2=C(NCCC3=CC=CC=C3)C(=O)C4=CC=CC=C4C2=O
CCC(=O)N(C(C)C)C1=C(NC)C(=O)C2=CC=CC=C2C1=O
CCN(CC)CCCC(C)NC1=CC=NC2=CC(Cl)=CC=C12
CC(CCCNCCO)NC1=CC=NC2=CC(Cl)=CC=C12
CC(C)C(=O)NCCCCNC1=CC=NC2=CC(Cl)=CC=C12
CC(C)(C)CC(=O)NCCCCNC1=CC=NC2=CC(Cl)=CC=C12
CC(C)CC(=O)NCCCCNC1=CC=NC2=CC(Cl)=CC=C12
CC(CC(=O)NCCCCNC1=CC=NC2=CC(Cl)=CC=C12)CC(C)(C)C
CCCCC(CC)C(=O)NCCCCNC1=CC=NC2=CC(Cl)=CC=C12
ClC1=CC=C2C(NCCCCNC(=O)C3CCCC3)=CC=NC2=C1
CCCCCCCCC(=O)NCCCCNC1=CC=NC2=CC(Cl)=CC=C12
CC(C)(CCl)C(=O)NCCCCNC1=CC=NC2=CC(Cl)=CC=C12
ClC1=CC=C2C(NCCCCNC(=O)C3CCCCC3)=CC=NC2=C1
NCCCCCCNC1=CC=NC2=CC(Cl)=CC=C12
ClC1=CC=C2C(NCCCCNC(=O)CCC3CCCC3)=CC=NC2=C1
ClC1=CC=C2C(NC3CCCCCCC3)=CC=NC2=C1
CC1CCC(CC1)NC2=CC=NC3=CC(Cl)=CC=C23
CN(C)C(=O)NCCCCNC1=CC=NC2=CC(Cl)=CC=C12
CCN(CC)CCCC(C)NC1=C2C=CC(Cl)=CC2=NC3=CC=C(OC)C=C13
CC(CCCO)NC1=CC=NC2=CC(Cl)=CC=C12
CCCCCC(=O)NCCCCCCNC1=CC=NC2=CC(Cl)=CC=C12
ClC1=CC=C2C(NC3CCC(CC3)NC(=O)CCC4CCCC4)=CC=NC2=C1
CN(C)CCCNC1=CC=NC2=CC(Cl)=CC=C12
'''

  mols = [Chem.MolFromSmiles(smiles) for smiles in molstr.splitlines()]
  na = len(mols)
  nb = len(mols)
  pathints = [getpathintegers(mol) for mol in mols]
  start = time.time()
  for a, api in zip(mols, pathints):
    for b, bpi in zip(mols, pathints):
      sim = AtomAtomPathSimilarity(a, b, m1pathintegers=api, m2pathintegers=bpi)
  print('time to compute %dx%d matrix: %.2fs' % (na, nb, time.time() - start))


class TestAtomAtomPathSimilarity(unittest.TestCase):

  def test_paper(self):
    self.assertEqual('%.3f' % test0(), '0.066')

  def test_getcommon(self):
    self.assertEqual(getcommon([2, 2, 2, 3, 3, 3], 6, [1, 2, 3, 3, 4, 5], 6), 3)

  def test_pathintegers(self):
    self.assertEqual(test1(), [{
      0: []
    }, {
      0: [1160, 2270],
      1: [2260, 30692],
      2: [1145, 33761]
    }, {
      0: [752, 1150, 1155, 3826, 38221, 43791],
      1: [1145, 1150, 1596, 4670, 32641, 38211],
      2: [1145, 1155, 5785, 32646, 43786, 65173]
    }])

  def test_AAPathComparator2Test(self):
    self.assertEqual(test2(), [
      '1.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000',
      '0.0000', '0.0000', '1.0000', '0.0000', '0.0345', '0.0182', '0.0000', '0.0000', '0.0000',
      '0.0017', '0.0000', '0.0000', '0.0000', '1.0000', '0.0000', '0.0182', '0.0000', '0.0000',
      '0.0000', '0.0020', '0.0000', '0.0000', '0.0345', '0.0000', '1.0000', '0.1126', '0.0000',
      '0.0000', '0.0000', '0.0088', '0.0000', '0.0000', '0.0182', '0.0182', '0.1126', '1.0000',
      '0.0000', '0.0000', '0.0000', '0.0336', '0.0000', '0.0000', '0.0000', '0.0000', '0.0000',
      '0.0000', '1.0000', '0.0373', '0.0645', '0.0148', '0.0869', '0.0000', '0.0000', '0.0000',
      '0.0000', '0.0000', '0.0373', '1.0000', '0.1101', '0.1767', '0.0869', '0.0000', '0.0000',
      '0.0000', '0.0000', '0.0000', '0.0645', '0.1101', '1.0000', '0.0387', '0.0219', '0.0000',
      '0.0017', '0.0020', '0.0088', '0.0336', '0.0148', '0.1767', '0.0387', '1.0000', '0.0869',
      '0.0000', '0.0000', '0.0000', '0.0000', '0.0000', '0.0869', '0.0869', '0.0219', '0.0869',
      '1.0000'
    ])

  def test_tableS1(self):
    self.assertEqual(test3(), [
      '1.00', '0.19', '0.06', '0.09', '0.19', '1.00', '0.12', '0.05', '0.06', '0.12', '1.00',
      '0.15', '0.09', '0.05', '0.15', '1.00'
    ])


if __name__ == "__main__":
  #	unittest.main()
  timeit()
