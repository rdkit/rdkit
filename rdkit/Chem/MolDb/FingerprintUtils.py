# $Id$
#
# Copyright (C) 2009 Greg Landrum
#  All Rights Reserved
#

import pickle

from rdkit import Chem, DataStructs

similarityMethods = {
  'RDK': DataStructs.ExplicitBitVect,
  'AtomPairs': DataStructs.IntSparseIntVect,
  'TopologicalTorsions': DataStructs.LongSparseIntVect,
  'Pharm2D': DataStructs.SparseBitVect,
  'Gobbi2D': DataStructs.SparseBitVect,
  'Morgan': DataStructs.UIntSparseIntVect,
  'Avalon': DataStructs.ExplicitBitVect,
}
supportedSimilarityMethods = list(iter(similarityMethods))


class LayeredOptions:
  loadLayerFlags = 0xFFFFFFFF
  searchLayerFlags = 0x7
  minPath = 1
  maxPath = 6
  fpSize = 1024
  wordSize = 32
  nWords = fpSize // wordSize

  @staticmethod
  def GetFingerprint(mol, query=True):
    if query:
      flags = LayeredOptions.searchLayerFlags
    else:
      flags = LayeredOptions.loadLayerFlags
    return Chem.LayeredFingerprint(mol, layerFlags=flags, minPath=LayeredOptions.minPath,
                                   maxPath=LayeredOptions.maxPath, fpSize=LayeredOptions.fpSize)

  @staticmethod
  def GetWords(mol, query=True):
    txt = LayeredOptions.GetFingerprint(mol, query=query).ToBitString()
    return [int(txt[x:x + 32], 2) for x in range(0, len(txt), 32)]

  @staticmethod
  def GetQueryText(mol, query=True):
    words = LayeredOptions.GetWords(mol, query=query)
    colqs = []
    for idx, word in enumerate(words):
      if not word:
        continue

      colqs.append(f'{word}&Col_{idx + 1}={word}')
    return ' and '.join(colqs)


def BuildSigFactory(options=None, fdefFile=None, bins=[(2, 3), (3, 4), (4, 5), (5, 6), (6, 7),
                                                       (7, 8), (8, 100)],
                    skipFeats=('LumpedHydrophobe', 'ZnBinder')):
  if options:
    fdefFile = options.fdefFile
  if not fdefFile:
    raise ValueError('bad fdef file')
  from rdkit.Chem import ChemicalFeatures
  from rdkit.Chem.Pharm2D import SigFactory
  featFactory = ChemicalFeatures.BuildFeatureFactory(fdefFile)
  sigFactory = SigFactory.SigFactory(featFactory, skipFeats=skipFeats, trianglePruneBins=False)
  sigFactory.SetBins(bins)
  return sigFactory


def BuildAtomPairFP(mol):
  from rdkit.Chem.AtomPairs import Pairs
  fp = Pairs.GetAtomPairFingerprintAsIntVect(mol)
  fp._sumCache = fp.GetTotalVal()
  return fp


def BuildTorsionsFP(mol):
  from rdkit.Chem.AtomPairs import Torsions
  fp = Torsions.GetTopologicalTorsionFingerprintAsIntVect(mol)
  fp._sumCache = fp.GetTotalVal()
  return fp


def BuildRDKitFP(mol):
  return Chem.RDKFingerprint(mol, nBitsPerHash=1)


def BuildPharm2DFP(mol):
  global sigFactory
  from rdkit.Chem.Pharm2D import Generate
  try:
    fp = Generate.Gen2DFingerprint(mol, sigFactory)
  except IndexError:
    print('FAIL:', Chem.MolToSmiles(mol, True))
    raise
  return fp


def BuildMorganFP(mol):
  from rdkit.Chem import rdMolDescriptors
  fp = rdMolDescriptors.GetMorganFingerprint(mol, 2)
  fp._sumCache = fp.GetTotalVal()
  return fp


def BuildAvalonFP(mol, smiles=None):
  from rdkit.Avalon import pyAvalonTools
  if smiles is None:
    return pyAvalonTools.GetAvalonFP(mol)
  return pyAvalonTools.GetAvalonFP(smiles, True)


def DepickleFP(pkl, similarityMethod):
  if not isinstance(pkl, (bytes, str)):
    pkl = str(pkl)
  try:
    klass = similarityMethods[similarityMethod]
    fp = klass(pkl)
  except Exception:
    import traceback
    traceback.print_exc()
    fp = pickle.loads(pkl)
  return fp
