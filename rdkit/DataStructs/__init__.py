# $Id$
#
#  Copyright (C) 2004-2006  Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#

from rdkit import rdBase
from rdkit.DataStructs import cDataStructs
from rdkit.DataStructs.cDataStructs import *

__doc__ = cDataStructs.__doc__

similarityFunctions = [
  ('Tanimoto', TanimotoSimilarity, ''),
  ("Dice", DiceSimilarity, ''),
  ("Cosine", CosineSimilarity, ''),
  ("Sokal", SokalSimilarity, ''),
  ("Russel", RusselSimilarity, ''),
  ("RogotGoldberg", RogotGoldbergSimilarity, ''),
  ("AllBit", AllBitSimilarity, ''),
  ("Kulczynski", KulczynskiSimilarity, ''),
  ("McConnaughey", McConnaugheySimilarity, ''),
  ("Asymmetric", AsymmetricSimilarity, ''),
  ("BraunBlanquet", BraunBlanquetSimilarity, ''),
]


def FingerprintSimilarity(fp1, fp2, metric=TanimotoSimilarity):
  """ returns the calculated similarity between two fingerprints,
      handles any folding that may need to be done to ensure that they
      are compatible

    """
  sz1 = fp1.GetNumBits()
  sz2 = fp2.GetNumBits()
  if sz1 < sz2:
    fp2 = FoldFingerprint(fp2, sz2 // sz1)
  elif sz2 < sz1:
    fp1 = FoldFingerprint(fp1, sz1 // sz2)
  return metric(fp1, fp2)


def FoldToTargetDensity(fp, density=0.3, minLength=64):
  while fp.GetNumOnBits() / len(fp) > density and len(fp) // 2 > minLength:
    fp = FoldFingerprint(fp, 2)
  return fp


ExplicitBitVect.ToBitString = BitVectToText
SparseBitVect.ToBitString = BitVectToText
