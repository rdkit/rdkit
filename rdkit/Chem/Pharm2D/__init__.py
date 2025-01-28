#
# Copyright (C) 2002-2009 greg Landrum and Rational Discovery LLC
#
""" module with functionality for 2D pharmacophores

"""


def DefaultSigFactory(fdefFile=None, minPointCount=2, maxPointCount=3, bins=[(2, 3), (3, 4), (4, 5),
                                                                             (5, 6), (6, 7), (7, 8),
                                                                             (8, 100)]):
  from rdkit.Chem import ChemicalFeatures
  from rdkit.Chem.Pharm2D import SigFactory
  if fdefFile is None:
    import os.path

    from rdkit import RDConfig
    fdefFile = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
  featFactory = ChemicalFeatures.BuildFeatureFactory(fdefFile, )
  factory = SigFactory.SigFactory(featFactory, skipFeats=('ZnBinder', 'LumpedHydrophobe'),
                                  minPointCount=minPointCount, maxPointCount=maxPointCount,
                                  trianglePruneBins=False)
  factory.SetBins(tuple(bins))
  factory.Init()
  return factory
