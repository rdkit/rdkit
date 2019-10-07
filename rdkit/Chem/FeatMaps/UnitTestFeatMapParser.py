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
from contextlib import closing
import unittest

from io import StringIO

from rdkit.Chem.FeatMaps import FeatMaps, FeatMapParser


def feq(n1, n2, tol=1e-5):
    return abs(n1 - n2) <= tol


class TestCase(unittest.TestCase):
    data = """

ScoreMode=Best
DirScoreMode=DotFullRange

BeginParams
  family=Aromatic radius=2.5 width=1.0 profile=Triangle
  family=Acceptor radius=1.5
EndParams

# optional
BeginPoints
  family=Acceptor pos=(1.0, 0.0, 5.0) weight=1.25 dir=(1, 1, 0)
  family=Aromatic pos=(0.0,1.0,0.0) weight=2.0 dir=(0,0,1) dir=(0,0,-1)
  family=Acceptor pos=(1.0,1.0,2.0) weight=1.25
EndPoints

"""

    def test1Basics(self):
        p = FeatMapParser.FeatMapParser()
        p.SetData(self.data)
        fm = p.Parse()
        self.assertTrue(fm.scoreMode == FeatMaps.FeatMapScoreMode.Best)
        self.assertTrue(fm.dirScoreMode == FeatMaps.FeatDirScoreMode.DotFullRange)
        self.assertTrue(fm.GetNumFeatures() == 3)

        feats = fm.GetFeatures()
        self.assertTrue(feq(feats[0].weight, 1.25))
        self.assertTrue(feq(feats[1].weight, 2.0))
        self.assertTrue(feq(feats[2].weight, 1.25))

        self.assertTrue(len(feats[0].featDirs) == 1)
        self.assertTrue(len(feats[1].featDirs) == 2)
        self.assertTrue(len(feats[2].featDirs) == 0)

        fams = [x.GetFamily() for x in feats]
        self.assertTrue(fams == ['Acceptor', 'Aromatic', 'Acceptor'])

    def test_FeatMapParser(self):
        # We can use a string
        p = FeatMapParser.FeatMapParser(data=self.data)
        fm = p.Parse()
        self.assertEqual(fm.GetNumFeatures(), 3)
        self.assertEqual([x.GetFamily() for x in fm.GetFeatures()],
                         ['Acceptor', 'Aromatic', 'Acceptor'])

        # We can use a list of strings
        p = FeatMapParser.FeatMapParser(data=self.data.split('\n'))
        fm = p.Parse()
        self.assertEqual(fm.GetNumFeatures(), 3)
        self.assertEqual([x.GetFamily() for x in fm.GetFeatures()],
                         ['Acceptor', 'Aromatic', 'Acceptor'])

        # and a stream
        with closing(StringIO(self.data)) as file:
            p = FeatMapParser.FeatMapParser(file=file)
        fm = p.Parse()
        self.assertEqual(fm.GetNumFeatures(), 3)
        self.assertEqual([x.GetFamily() for x in fm.GetFeatures()],
                         ['Acceptor', 'Aromatic', 'Acceptor'])

    def test_ParseErrors(self):
        # Typos in scoreMode or dirscoreMode section
        data = "scoreMode = typo\nbeginParams\nfamily=Acceptor radius=1.5\nEndParams"
        p = FeatMapParser.FeatMapParser(data=data)
        self.assertRaises(FeatMapParser.FeatMapParseError, p.Parse)

        data = "dirscoremode = typo\nbeginParams\nfamily=Acceptor radius=1.5\nEndParams"
        p = FeatMapParser.FeatMapParser(data=data)
        self.assertRaises(FeatMapParser.FeatMapParseError, p.Parse)

        data = "typo = All\nbeginParams\nfamily=Acceptor radius=1.5\nEndParams"
        p = FeatMapParser.FeatMapParser(data=data)
        self.assertRaises(FeatMapParser.FeatMapParseError, p.Parse)

        # Typos in paramBlock
        data = "beginTypo\nfamily=Acceptor radius=1.5\nEndParams"
        p = FeatMapParser.FeatMapParser(data=data)
        self.assertRaises(FeatMapParser.FeatMapParseError, p.Parse)

        data = "beginParams\nfamily=Acceptor radius=1.5\nEndTypo"
        p = FeatMapParser.FeatMapParser(data=data)
        self.assertRaises(FeatMapParser.FeatMapParseError, p.Parse)

        data = "beginParams\ntypo=Acceptor radius=1.5\nEndParams"
        p = FeatMapParser.FeatMapParser(data=data)
        self.assertRaises(FeatMapParser.FeatMapParseError, p.Parse)

        data = "beginParams\nprofile=Typo\nEndParams"
        p = FeatMapParser.FeatMapParser(data=data)
        self.assertRaises(FeatMapParser.FeatMapParseError, p.Parse)

        # Typos in points block
        data = "BeginPoints\npos=(1.0, 0.0, 5.0, 4.0)\nEndPoints"
        p = FeatMapParser.FeatMapParser(data=data)
        self.assertRaises(ValueError, p.Parse)

        data = "BeginPoints\npos=(1.0, 0.0, 5.0) typo=Acceptor\nEndPoints"
        p = FeatMapParser.FeatMapParser(data=data)
        self.assertRaises(FeatMapParser.FeatMapParseError, p.Parse)


if __name__ == '__main__':  # pragma: nocover
    unittest.main()
