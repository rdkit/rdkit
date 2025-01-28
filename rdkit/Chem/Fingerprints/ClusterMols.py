#
#  Copyright (c) 2003-2006 Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" utility functionality for clustering molecules using fingerprints
 includes a command line app for clustering


Sample Usage:
  python ClusterMols.py  -d data.gdb -t daylight_sig \
    --idName="CAS_TF" -o clust1.pkl \
    --actTable="dop_test" --actName="moa_quant"

"""

import pickle

import numpy

from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols, MolSimilarity
from rdkit.ML.Cluster import Murtagh

message = FingerprintMols.message
error = FingerprintMols.error


def GetDistanceMatrix(data, metric, isSimilarity=1):
  """ data should be a list of tuples with fingerprints in position 1
   (the rest of the elements of the tuple are not important)

    Returns the symmetric distance matrix
    (see ML.Cluster.Resemblance for layout documentation)

  """
  nPts = len(data)
  distsMatrix = numpy.zeros((nPts * (nPts - 1) // 2), dtype=numpy.float64)
  nSoFar = 0
  for col in range(1, nPts):
    fp1 = data[col][1]
    nBits1 = fp1.GetNumBits()
    for row in range(col):
      fp2 = data[row][1]
      nBits2 = fp2.GetNumBits()
      if nBits1 > nBits2:
        fp1 = DataStructs.FoldFingerprint(fp1, nBits1 / nBits2)
      elif nBits2 > nBits1:
        fp2 = DataStructs.FoldFingerprint(fp2, nBits2 / nBits1)

      if isSimilarity:
        distsMatrix[nSoFar] = 1.0 - metric(fp1, fp2)
      else:
        distsMatrix[nSoFar] = metric(fp1, fp2)
      nSoFar += 1
  return distsMatrix


def ClusterPoints(data, metric, algorithmId, haveLabels=False, haveActs=True,
                  returnDistances=False):
  message('Generating distance matrix.\n')
  dMat = GetDistanceMatrix(data, metric)
  message('Clustering\n')
  clustTree = Murtagh.ClusterData(dMat, len(data), algorithmId, isDistData=1)[0]
  acts = []
  if haveActs and len(data[0]) > 2:
    # we've got activities... use them:
    acts = [int(x[2]) for x in data]

  if not haveLabels:
    labels = [f'Mol: {x[0]}' for x in data]
  else:
    labels = [x[0] for x in data]
  clustTree._ptLabels = labels
  if acts:
    clustTree._ptValues = acts

  for pt in clustTree.GetPoints():
    idx = pt.GetIndex() - 1
    pt.SetName(labels[idx])
    if acts:
      try:
        pt.SetData(int(acts[idx]))
      except Exception:
        pass

  if not returnDistances:
    return clustTree
  return clustTree, dMat


def ClusterFromDetails(details):
  """ Returns the cluster tree

  """
  data = MolSimilarity.GetFingerprints(details)
  if details.maxMols > 0:
    data = data[:details.maxMols]
  if details.outFileName:
    try:
      outF = open(details.outFileName, 'wb+')
    except IOError:
      error("Error: could not open output file %s for writing\n" % (details.outFileName))
      return None
  else:
    outF = None

  if not data:
    return None

  clustTree = ClusterPoints(data, details.metric, details.clusterAlgo, haveLabels=0, haveActs=1)
  if outF:
    pickle.dump(clustTree, outF)
  return clustTree


_usageDoc = """
Usage: ClusterMols.py [args] <fName>

  If <fName> is provided and no tableName is specified (see below),
  data will be read from the text file <fName>.  Text files delimited
  with either commas (extension .csv) or tabs (extension .txt) are
  supported.

  Command line arguments are:

    - -d _dbName_: set the name of the database from which
      to pull input fingerprint information.

    - -t _tableName_: set the name of the database table
      from which to pull input fingerprint information

    - --idName=val: sets the name of the id column in the input
      database.  Default is *ID*.

    - -o _outFileName_:  name of the output file (output will
      be a pickle (.pkl) file with the cluster tree)

    - --actTable=val: name of table containing activity values
     (used to color points in the cluster tree).

    - --actName=val: name of column with activities in the activity
      table.  The values in this column should either be integers or
      convertible into integers.

    - --SLINK: use the single-linkage clustering algorithm
      (default is Ward's minimum variance)

    - --CLINK: use the complete-linkage clustering algorithm
      (default is Ward's minimum variance)

    - --UPGMA: use the group-average clustering algorithm
      (default is Ward's minimum variance)

    - --dice: use the DICE similarity metric instead of Tanimoto

    - --cosine: use the cosine similarity metric instead of Tanimoto

    - --fpColName=val: name to use for the column which stores
      fingerprints (in pickled format) in the input db table.
      Default is *AutoFragmentFP*

    - --minPath=val:  minimum path length to be included in
      fragment-based fingerprints. Default is *2*.

    - --maxPath=val:  maximum path length to be included in
      fragment-based fingerprints. Default is *7*.

    - --nBitsPerHash: number of bits to be set in the output
      fingerprint for each fragment. Default is *4*.

    - --discrim: use of path-based discriminators to hash bits.
      Default is *false*.

    - -V: include valence information in the fingerprints
      Default is *false*.

    - -H: include Hs in the fingerprint
      Default is *false*.

    - --useMACCS: use the public MACCS keys to do the fingerprinting
      (instead of a daylight-type fingerprint)


"""
if __name__ == '__main__':
  message("This is ClusterMols\n\n")
  FingerprintMols._usageDoc = _usageDoc
  details = FingerprintMols.ParseArgs()
  ClusterFromDetails(details)
