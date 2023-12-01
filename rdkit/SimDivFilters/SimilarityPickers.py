# $Id$
#
# Copyright (C) 2003-2008 Greg Landrum and Rational Discovery LLC
#  All Rights Reserved
#

import bisect

from rdkit import DataStructs
from rdkit.DataStructs.TopNContainer import TopNContainer


class GenericPicker(object):
  _picks = None

  def MakePicks(self, force=False):
    raise NotImplementedError("GenericPicker is a virtual base class")

  def __len__(self):
    if self._picks is None:
      self.MakePicks()
    return len(self._picks)

  def __getitem__(self, which):
    if self._picks is None:
      self.MakePicks()
    return self._picks[which]


class TopNOverallPicker(GenericPicker):
  """  A class for picking the top N overall best matches across a library

  Connect to a database and build molecules:

  >>> from rdkit import Chem
  >>> from rdkit import RDConfig
  >>> import os.path
  >>> from rdkit.Dbase.DbConnection import DbConnect
  >>> dbName = RDConfig.RDTestDatabase
  >>> conn = DbConnect(dbName,'simple_mols1')
  >>> [x.upper() for x in conn.GetColumnNames()]
  ['SMILES', 'ID']
  >>> mols = []
  >>> for smi,id in conn.GetData():
  ...   mol = Chem.MolFromSmiles(str(smi))
  ...   mol.SetProp('_Name',str(id))
  ...   mols.append(mol)
  >>> len(mols)
  12

  Calculate fingerprints:

  >>> probefps = []
  >>> for mol in mols:
  ...   fp = Chem.RDKFingerprint(mol)
  ...   fp._id = mol.GetProp('_Name')
  ...   probefps.append(fp)

  Start by finding the top matches for a single probe.  This ether should pull
  other ethers from the db:

  >>> mol = Chem.MolFromSmiles('COC')
  >>> probeFp = Chem.RDKFingerprint(mol)
  >>> picker = TopNOverallPicker(numToPick=2,probeFps=[probeFp],dataSet=probefps)
  >>> len(picker)
  2
  >>> fp,score = picker[0]
  >>> id = fp._id
  >>> str(id)
  'ether-1'
  >>> score
  1.0

  The results come back in order:

  >>> fp,score = picker[1]
  >>> id = fp._id
  >>> str(id)
  'ether-2'

  Now find the top matches for 2 probes.  We'll get one ether and one acid:

  >>> fps = []
  >>> fps.append(Chem.RDKFingerprint(Chem.MolFromSmiles('COC')))
  >>> fps.append(Chem.RDKFingerprint(Chem.MolFromSmiles('CC(=O)O')))
  >>> picker = TopNOverallPicker(numToPick=3,probeFps=fps,dataSet=probefps)
  >>> len(picker)
  3
  >>> fp,score = picker[0]
  >>> id = fp._id
  >>> str(id)
  'acid-1'
  >>> fp,score = picker[1]
  >>> id = fp._id
  >>> str(id)
  'ether-1'
  >>> score
  1.0
  >>> fp,score = picker[2]
  >>> id = fp._id
  >>> str(id)
  'acid-2'

  """

  def __init__(self, numToPick=10, probeFps=None, dataSet=None,
               simMetric=DataStructs.TanimotoSimilarity):
    """

      dataSet should be a sequence of BitVectors

    """
    self.numToPick = numToPick
    self.probes = probeFps
    self.data = dataSet
    self.simMetric = simMetric
    self._picks = None

  def MakePicks(self, force=False):
    if self._picks is not None and not force:
      return
    picks = TopNContainer(self.numToPick)
    for fp in self.data:
      origFp = fp
      bestScore = -1.0
      for probeFp in self.probes:
        score = DataStructs.FingerprintSimilarity(origFp, probeFp, self.simMetric)
        bestScore = max(score, bestScore)
      picks.Insert(bestScore, fp)
    self._picks = []
    for score, pt in picks:
      self._picks.append((pt, score))
    self._picks.reverse()


class SpreadPicker(GenericPicker):
  """  A class for picking the best matches across a library

  Connect to a database:

  >>> from rdkit import Chem
  >>> from rdkit import RDConfig
  >>> import os.path
  >>> from rdkit.Dbase.DbConnection import DbConnect
  >>> dbName = RDConfig.RDTestDatabase
  >>> conn = DbConnect(dbName,'simple_mols1')
  >>> [x.upper() for x in conn.GetColumnNames()]
  ['SMILES', 'ID']
  >>> mols = []
  >>> for smi,id in conn.GetData():
  ...   mol = Chem.MolFromSmiles(str(smi))
  ...   mol.SetProp('_Name',str(id))
  ...   mols.append(mol)
  >>> len(mols)
  12

  Calculate fingerprints:

  >>> probefps = []
  >>> for mol in mols:
  ...   fp = Chem.RDKFingerprint(mol)
  ...   fp._id = mol.GetProp('_Name')
  ...   probefps.append(fp)

  Start by finding the top matches for a single probe.  This ether should pull
  other ethers from the db:

  >>> mol = Chem.MolFromSmiles('COC')
  >>> probeFp = Chem.RDKFingerprint(mol)
  >>> picker = SpreadPicker(numToPick=2,probeFps=[probeFp],dataSet=probefps)
  >>> len(picker)
  2
  >>> fp,score = picker[0]
  >>> id = fp._id
  >>> str(id)
  'ether-1'
  >>> score
  1.0

  The results come back in order:

  >>> fp,score = picker[1]
  >>> id = fp._id
  >>> str(id)
  'ether-2'

  Now find the top matches for 2 probes.  We'll get one ether and one acid:

  >>> fps = []
  >>> fps.append(Chem.RDKFingerprint(Chem.MolFromSmiles('COC')))
  >>> fps.append(Chem.RDKFingerprint(Chem.MolFromSmiles('CC(=O)O')))
  >>> picker = SpreadPicker(numToPick=3,probeFps=fps,dataSet=probefps)
  >>> len(picker)
  3
  >>> fp,score = picker[0]
  >>> id = fp._id
  >>> str(id)
  'ether-1'
  >>> score
  1.0
  >>> fp,score = picker[1]
  >>> id = fp._id
  >>> str(id)
  'acid-1'
  >>> score
  1.0
  >>> fp,score = picker[2]
  >>> id = fp._id
  >>> str(id)
  'ether-2'

  """

  def __init__(self, numToPick=10, probeFps=None, dataSet=None,
               simMetric=DataStructs.TanimotoSimilarity, expectPickles=True, onlyNames=False):
    """

      dataSet should be a sequence of BitVectors or, if expectPickles
      is False, a set of strings that can be converted to bit vectors

    """
    self.numToPick = numToPick
    self.probes = probeFps
    self.data = dataSet
    self.simMetric = simMetric
    self.expectPickles = expectPickles
    self.onlyNames = onlyNames

    self._picks = None

  def MakePicks(self, force=False, silent=False):
    if self._picks is not None and not force:
      return

    # start by getting the NxM score matrix
    #  (N=num probes, M=num fps)
    nProbes = len(self.probes)
    scores = [None] * nProbes
    for i in range(nProbes):
      scores[i] = []
    j = 0
    fps = []
    for origFp in self.data:
      for i in range(nProbes):
        score = DataStructs.FingerprintSimilarity(self.probes[i], origFp, self.simMetric)
        bisect.insort(scores[i], (score, j))
        if len(scores[i]) >= self.numToPick:
          del scores[self.numToPick:]
      if self.onlyNames and hasattr(origFp, '_fieldsFromDb'):
        fps.append(origFp._fieldsFromDb[0])
      else:
        fps.append(origFp)
      j += 1
      if not silent and not j % 1000:
        print('scored %d fps' % j)

    # now go probe by probe and select the current top entry until we are finished:
    nPicked = 0
    self._picks = []
    taken = [0] * len(fps)
    while nPicked < self.numToPick:
      rowIdx = nPicked % len(scores)
      row = scores[rowIdx]
      score, idx = row.pop()
      # make sure we haven't taken this one already (from another row):
      while taken[idx] and len(row):
        score, idx = row.pop()
      if not taken[idx]:
        fp = fps[idx]
        self._picks.append((fp, score))
        taken[idx] = 1
        nPicked += 1


# ------------------------------------
#
#  doctest boilerplate
#
def _runDoctests(verbose=None):  # pragma: nocover
  import doctest
  import sys
  failed, _ = doctest.testmod(optionflags=doctest.ELLIPSIS, verbose=verbose)
  sys.exit(failed)


if __name__ == '__main__':  # pragma: nocover
  _runDoctests()
