# $Id$
#
#  Copyright (C) 2006  greg Landrum
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
import doctest
import math
import unittest

from rdkit.Chem.ChemicalFeatures import FreeChemicalFeature
from rdkit.Chem.FeatMaps import FeatMaps
from rdkit.Geometry import Point3D


def feq(n1, n2, tol=1e-5):
  return abs(n1 - n2) <= tol


def load_tests(loader, tests, ignore):
  """ Add the Doctests from the module """
  tests.addTests(doctest.DocTestSuite(FeatMaps, optionflags=doctest.ELLIPSIS))
  return tests


class TestCase(unittest.TestCase):

  def test1Basics(self):
    fs = [FreeChemicalFeature('Aromatic', 'Foo', Point3D(0, 0, 0)),
          FreeChemicalFeature('Acceptor', 'Foo', Point3D(2, 0, 0)),
          FreeChemicalFeature('Acceptor', 'Foo', Point3D(2.1, 0, 0))]
    aFmp = FeatMaps.FeatMapParams()
    aFmp.radius = 2.0
    bFmp = FeatMaps.FeatMapParams()
    bFmp.radius = 1.0

    fmps = {'Aromatic': aFmp, 'Acceptor': bFmp}

    fmap = FeatMaps.FeatMap(params=fmps)
    fmap.AddFeature(fs[0], 1.0)
    fmap.AddFeature(fs[1], 1.0)
    fmap.AddFeature(fs[2], 2.0)
    fmap.AddFeature(fs[2], 2.0)

    self.assertTrue(fmap.GetNumFeatures() == 4)
    self.assertTrue(len(fmap.GetFeatures()) == 4)
    fmap.DropFeature(3)
    self.assertTrue(fmap.GetNumFeatures() == 3)
    self.assertTrue(len(fmap.GetFeatures()) == 3)

    f = fmap.GetFeature(0)
    self.assertTrue(f.GetFamily() == 'Aromatic')
    self.assertTrue(feq(f.weight, 1.0))
    f = fmap.GetFeature(1)
    self.assertTrue(f.GetFamily() == 'Acceptor')
    self.assertTrue(feq(f.weight, 1.0))
    f = fmap.GetFeature(2)
    self.assertTrue(f.GetFamily() == 'Acceptor')
    self.assertTrue(feq(f.weight, 2.0))

  def test2FeatFeatScoreGauss(self):
    aFmp = FeatMaps.FeatMapParams()
    aFmp.radius = 2.0
    fmps = {'Aromatic': aFmp}
    fmap = FeatMaps.FeatMap(params=fmps)

    fs = [FreeChemicalFeature('Aromatic', 'Foo', Point3D(0, 0, 0))]
    fmap.AddFeature(fs[0], 1.0)
    self.assertTrue(len(fmap.GetFeatures()) == 1)

    sc = fmap.GetFeatFeatScore(
      fmap.GetFeature(0), FreeChemicalFeature('Aromatic', '', Point3D(1, 0, 0)))
    self.assertTrue(feq(sc, math.exp(-1)))
    sc = fmap.GetFeatFeatScore(
      fmap.GetFeature(0), FreeChemicalFeature('Aromatic', '', Point3D(1.5, 0, 0)))
    self.assertTrue(feq(sc, math.exp(-2.25)))
    sc = fmap.GetFeatFeatScore(
      fmap.GetFeature(0), FreeChemicalFeature('Aromatic', '', Point3D(0, 0, 0)))
    self.assertTrue(feq(sc, 1.0))
    sc = fmap.GetFeatFeatScore(
      fmap.GetFeature(0), FreeChemicalFeature('Aromatic', '', Point3D(2.1, 0, 0)))
    self.assertTrue(feq(sc, 0))

    sc = fmap.GetFeatFeatScore(
      fmap.GetFeature(0), FreeChemicalFeature('Acceptor', '', Point3D(1, 0, 0)))
    self.assertTrue(feq(sc, 0))
    sc = fmap.GetFeatFeatScore(
      fmap.GetFeature(0), FreeChemicalFeature('Acceptor', '', Point3D(1, 0, 0)), typeMatch=False)
    self.assertTrue(feq(sc, math.exp(-1)))

    self.assertRaises(IndexError,lambda: fmap.GetFeatFeatScore(fmap.GetFeature(1),FreeChemicalFeature('Aromatic','',Point3D(0,0,0))))

  def test3FeatFeatScoreTriangle(self):
    aFmp = FeatMaps.FeatMapParams()
    aFmp.width = 2.0
    aFmp.radius = 3.0
    aFmp.featProfile = FeatMaps.FeatMapParams.FeatProfile.Triangle
    fmps = {'Aromatic': aFmp}
    fmap = FeatMaps.FeatMap(params=fmps)

    fs = [FreeChemicalFeature('Aromatic', 'Foo', Point3D(0, 0, 0))]
    fmap.AddFeature(fs[0], 1.0)
    self.assertTrue(len(fmap.GetFeatures()) == 1)

    sc = fmap.GetFeatFeatScore(
      fmap.GetFeature(0), FreeChemicalFeature('Aromatic', '', Point3D(1, 0, 0)))
    self.assertTrue(feq(sc, 0.5))
    sc = fmap.GetFeatFeatScore(
      fmap.GetFeature(0), FreeChemicalFeature('Aromatic', '', Point3D(1.5, 0, 0)))
    self.assertTrue(feq(sc, 0.25))
    sc = fmap.GetFeatFeatScore(
      fmap.GetFeature(0), FreeChemicalFeature('Aromatic', '', Point3D(0, 0, 0)))
    self.assertTrue(feq(sc, 1.0))
    sc = fmap.GetFeatFeatScore(
      fmap.GetFeature(0), FreeChemicalFeature('Aromatic', '', Point3D(2.1, 0, 0)))
    self.assertTrue(feq(sc, 0))

  def test4FeatFeatScoreBox(self):
    aFmp = FeatMaps.FeatMapParams()
    aFmp.radius = 2.0
    aFmp.featProfile = FeatMaps.FeatMapParams.FeatProfile.Box
    fmps = {'Aromatic': aFmp}
    fmap = FeatMaps.FeatMap(params=fmps)

    fs = [FreeChemicalFeature('Aromatic', 'Foo', Point3D(0, 0, 0))]
    fmap.AddFeature(fs[0], 1.1)
    self.assertTrue(len(fmap.GetFeatures()) == 1)

    sc = fmap.GetFeatFeatScore(
      fmap.GetFeature(0), FreeChemicalFeature('Aromatic', '', Point3D(1, 0, 0)))
    self.assertTrue(feq(sc, 1.1))
    sc = fmap.GetFeatFeatScore(
      fmap.GetFeature(0), FreeChemicalFeature('Aromatic', '', Point3D(1.5, 0, 0)))
    self.assertTrue(feq(sc, 1.1))
    sc = fmap.GetFeatFeatScore(
      fmap.GetFeature(0), FreeChemicalFeature('Aromatic', '', Point3D(0, 0, 0)))
    self.assertTrue(feq(sc, 1.1))
    sc = fmap.GetFeatFeatScore(
      fmap.GetFeature(0), FreeChemicalFeature('Aromatic', '', Point3D(2.1, 0, 0)))
    self.assertTrue(feq(sc, 0))

  def test5ScoreFeats(self):
    aFmp = FeatMaps.FeatMapParams()
    bFmp = FeatMaps.FeatMapParams()
    fmps = {'Aromatic': aFmp, 'Acceptor': bFmp}
    fmap = FeatMaps.FeatMap(params=fmps)

    fs = [FreeChemicalFeature('Aromatic', 'Foo', Point3D(0, 0, 0)),
          FreeChemicalFeature('Acceptor', 'Foo', Point3D(2, 0, 0)),
          FreeChemicalFeature('Acceptor', 'Foo', Point3D(2.1, 0, 0))]
    fmap.AddFeature(fs[0], 1.1)
    fmap.AddFeature(fs[1], 1.1)
    fmap.AddFeature(fs[2], 2.1)

    l1 = fmap._loopOverMatchingFeats(FreeChemicalFeature('Aromatic', '', Point3D(0, 0, 0)))
    l1 = list(l1)
    self.assertTrue(len(l1) == 1)
    self.assertTrue(l1[0][0] == 0)
    self.assertTrue(l1[0][1].GetFamily() == 'Aromatic')

    l1 = fmap._loopOverMatchingFeats(FreeChemicalFeature('Acceptor', '', Point3D(0, 0, 0)))
    l1 = list(l1)
    self.assertTrue(len(l1) == 2)
    self.assertTrue(l1[0][0] == 1)
    self.assertTrue(l1[0][1].GetFamily() == 'Acceptor')
    self.assertTrue(l1[1][0] == 2)
    self.assertTrue(l1[1][1].GetFamily() == 'Acceptor')

  def test6ScoreFeats(self):
    aFmp = FeatMaps.FeatMapParams()
    aFmp.radius = 2.0
    bFmp = FeatMaps.FeatMapParams()
    bFmp.radius = 1.0

    # use box scoring to make math easier... the other scoring
    # functions are tested above
    aFmp.featProfile = FeatMaps.FeatMapParams.FeatProfile.Box
    bFmp.featProfile = FeatMaps.FeatMapParams.FeatProfile.Box

    fmps = {'Aromatic': aFmp, 'Acceptor': bFmp}

    fmap = FeatMaps.FeatMap(params=fmps)

    fs = [FreeChemicalFeature('Aromatic', 'Foo', Point3D(0, 0, 0)),
          FreeChemicalFeature('Acceptor', 'Foo', Point3D(2, 0, 0)),
          FreeChemicalFeature('Acceptor', 'Foo', Point3D(2.1, 0, 0))]
    fmap.AddFeature(fs[0], 1.1)
    fmap.AddFeature(fs[1], 1.1)
    fmap.AddFeature(fs[2], 2.1)

    fs = [FreeChemicalFeature('Aromatic', '', Point3D(0, 1, 0)),
          FreeChemicalFeature('Acceptor', '', Point3D(1.5, 0, 0)), ]

    sc = fmap.ScoreFeats(fs)
    self.assertTrue(feq(sc, 4.3))

    msv = [-1] * 3
    fsv = [-1] * 2
    fsfmi = [None] * 2
    sc = fmap.ScoreFeats(fs, mapScoreVect=msv, featsScoreVect=fsv, featsToFeatMapIdx=fsfmi)
    self.assertTrue(feq(sc, 4.3))
    self.assertTrue(feq(sum(msv), sc))
    self.assertTrue(feq(sum(fsv), sc))
    self.assertTrue(fsfmi == [[0], [1, 2]])

    # make sure we reset the vectors internally:
    sc = fmap.ScoreFeats(fs, mapScoreVect=msv, featsScoreVect=fsv, featsToFeatMapIdx=fsfmi)
    self.assertTrue(feq(sc, 4.3))
    self.assertTrue(feq(sum(msv), sc))
    self.assertTrue(feq(sum(fsv), sc))
    self.assertTrue(fsfmi == [[0], [1, 2]])

    fmap.scoreMode = FeatMaps.FeatMapScoreMode.Closest
    sc = fmap.ScoreFeats(fs, mapScoreVect=msv, featsScoreVect=fsv, featsToFeatMapIdx=fsfmi)
    self.assertTrue(feq(sc, 2.1))
    self.assertTrue(feq(sum(msv), sc))
    self.assertTrue(feq(sum(fsv), sc))
    self.assertTrue(fsfmi == [[0], [1]])

    fmap.scoreMode = FeatMaps.FeatMapScoreMode.Best
    sc = fmap.ScoreFeats(fs, mapScoreVect=msv, featsScoreVect=fsv, featsToFeatMapIdx=fsfmi)
    self.assertTrue(feq(sc, 3.2))
    self.assertTrue(feq(sum(msv), sc))
    self.assertTrue(feq(sum(fsv), sc))
    self.assertTrue(fsfmi == [[0], [2]])

  def test7ScoreFeats(self):
    aFmp = FeatMaps.FeatMapParams()
    aFmp.radius = 2.0
    aFmp.featProfile = FeatMaps.FeatMapParams.FeatProfile.Box
    bFmp = FeatMaps.FeatMapParams()
    bFmp.radius = 1.0
    bFmp.featProfile = FeatMaps.FeatMapParams.FeatProfile.Box

    fmps = {'Aromatic': aFmp, 'Acceptor': bFmp}
    fmap = FeatMaps.FeatMap(params=fmps)

    fs = [FreeChemicalFeature('Aromatic', 'Foo', Point3D(0, 0, 0)), ]
    fmap.AddFeature(fs[0], 1.1)

    fs = [FreeChemicalFeature('Aromatic', '', Point3D(0, 1, 0)),
          FreeChemicalFeature('Acceptor', '', Point3D(1.5, 0, 0)), ]

    sc = fmap.ScoreFeats(fs)
    self.assertTrue(feq(sc, 1.1))

    msv = [-1] * 1
    fsv = [-1] * 2
    fsfmi = [-1] * 2
    sc = fmap.ScoreFeats(fs, mapScoreVect=msv, featsScoreVect=fsv, featsToFeatMapIdx=fsfmi)
    self.assertTrue(feq(sc, 1.1))
    self.assertTrue(feq(sum(msv), sc))
    self.assertTrue(feq(sum(fsv), sc))
    self.assertTrue(fsfmi == [[0], []])

    fmap.scoreMode = FeatMaps.FeatMapScoreMode.Closest
    sc = fmap.ScoreFeats(fs, mapScoreVect=msv, featsScoreVect=fsv, featsToFeatMapIdx=fsfmi)
    self.assertTrue(feq(sc, 1.1))
    self.assertTrue(feq(sum(msv), sc))
    self.assertTrue(feq(sum(fsv), sc))
    self.assertTrue(fsfmi == [[0], []])

    fmap.scoreMode = FeatMaps.FeatMapScoreMode.Best
    sc = fmap.ScoreFeats(fs, mapScoreVect=msv, featsScoreVect=fsv, featsToFeatMapIdx=fsfmi)
    self.assertTrue(feq(sc, 1.1))
    self.assertTrue(feq(sum(msv), sc))
    self.assertTrue(feq(sum(fsv), sc))
    self.assertTrue(fsfmi == [[0], []])

  def test8ScoreFeatDirs(self):
    aFmp = FeatMaps.FeatMapParams()
    aFmp.radius = 2.0
    aFmp.featProfile = FeatMaps.FeatMapParams.FeatProfile.Box

    fmps = {'Acceptor': aFmp}
    fmap = FeatMaps.FeatMap(params=fmps)
    fmap.dirScoreMode = FeatMaps.FeatDirScoreMode.DotFullRange

    fs = [FreeChemicalFeature('Acceptor', 'Foo', Point3D(0, 0, 0)), ]
    fs[0].featDirs = [Point3D(0, 1, 0)]
    fmap.AddFeature(fs[0], 1.1)

    fs = [FreeChemicalFeature('Acceptor', '', Point3D(1.0, 0, 0)), ]
    fs[0].featDirs = [Point3D(0, 1, 0)]

    sc = fmap.ScoreFeats(fs)
    self.assertTrue(feq(sc, 1.1))

    fs[0].featDirs = [Point3D(0, -1, 0)]
    sc = fmap.ScoreFeats(fs)
    self.assertTrue(feq(sc, -1.1))

    fs[0].featDirs = [Point3D(0, 0, 1)]
    sc = fmap.ScoreFeats(fs)
    self.assertTrue(feq(sc, 0.0))

    fmap.dirScoreMode = FeatMaps.FeatDirScoreMode.DotPosRange
    fs[0].featDirs = [Point3D(0, -1, 0)]
    sc = fmap.ScoreFeats(fs)
    self.assertTrue(feq(sc, 0))

    fs[0].featDirs = [Point3D(0, 1, 0)]
    sc = fmap.ScoreFeats(fs)
    self.assertTrue(feq(sc, 1.1))


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
