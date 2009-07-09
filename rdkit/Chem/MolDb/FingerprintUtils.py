# $Id$
#
# Copyright (C) 2009 Greg Landrum
#  All Rights Reserved
#
import cPickle
from rdkit import DataStructs

def BuildSigFactory(options=None,fdefFile=None,
                    bins=[(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,100)],
                    skipFeats=('LumpedHydrophobe','ZnBinder')):
  if options:
    fdefFile = options.fdefFile
  if not fdefFile:
    raise ValueError,'bad fdef file'
  from rdkit.Chem import ChemicalFeatures
  from rdkit.Chem.Pharm2D import SigFactory
  featFactory = ChemicalFeatures.BuildFeatureFactory(fdefFile)
  sigFactory = SigFactory.SigFactory(featFactory,
                                     skipFeats=skipFeats,
                                     trianglePruneBins=False)
  sigFactory.SetBins(bins)
  return sigFactory

similarityMethods={'RDK':DataStructs.ExplicitBitVect,
                   'AtomPairs':DataStructs.IntSparseIntVect,
                   'TopologicalTorsions':DataStructs.LongSparseIntVect,
                   'Pharm2D':DataStructs.SparseBitVect,
                   'Gobbi2D':DataStructs.SparseBitVect,
                   'Morgan':DataStructs.UIntSparseIntVect
                   }
supportedSimilarityMethods=similarityMethods.keys()


def DepickleFP(pkl,similarityMethod):
    try:
        klass=similarityMethods[similarityMethod]
        fp = klass(pkl)
    except:
        import traceback
        traceback.print_exc()
        fp = cPickle.loads(pkl)
    return fp
