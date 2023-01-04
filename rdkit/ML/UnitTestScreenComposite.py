# $Id$
#
#  Copyright (C) 2003-2008  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""unit testing code for the ScreenComposite functionality

"""
import io
import os
import unittest

from rdkit import RDConfig
from rdkit.ML import ScreenComposite
import pickle
VERSION = 2
DUMP = False

class TestCase(unittest.TestCase):

    def setUp(self):
        self.baseDir = os.path.join(RDConfig.RDCodeDir, 'ML', 'test_data')
        self.baseDir = os.path.join(os.path.dirname(__file__), 'test_data')
        self.dbName = RDConfig.RDTestDatabase
        self.details = ScreenComposite.SetDefaults()
        self.details.dbName = self.dbName
        self.details.dbUser = RDConfig.defaultDBUser
        self.details.dbPassword = RDConfig.defaultDBPassword

    def _open(self, path):
        if VERSION == 1:
            with open(path, 'r') as pklTF:
                buf = pklTF.read().replace('\r\n', '\n').encode('utf-8')
                pklTF.close()
            with io.BytesIO(buf) as pklF:
                return pickle.load(pklF)
        else:
            with open(path, 'rb') as buf:
                return pickle.load(buf)

    def _dump(self,n, nGood, misCount, nSkipped, avgGood, avgBad, avgSkip, tbl):
        if not DUMP: return
        print(f"""
        test {n}
        self.assertEqual(nGood, {nGood})
        self.assertEqual(misCount, {misCount})
        self.assertEqual(nSkipped, {nSkipped})
        self.assertAlmostEqual(avgGood, {avgGood:0.4f}, 4)
        self.assertAlmostEqual(avgBad, {avgBad:0.4f}, 4)
        self.assertAlmostEqual(avgSkip, {avgSkip:0.4f}, 4)
        self.assertEqual(tbl[0, 0], {tbl[0,0]})
        self.assertEqual(tbl[1, 1], {tbl[1,1]})
        self.assertEqual(tbl[0, 1], {tbl[0,1]})
        self.assertEqual(tbl[1, 0], {tbl[1,0]})
        """)
        
    def test1_basics(self):
        # """ basics """
        self.details.tableName = 'ferro_quant'
        compos = self._open(os.path.join(self.baseDir, 'ferromag_quant_10.pkl'))
        tgt = 7
        self.assertEqual(len(compos), tgt)

        nGood, misCount, nSkipped, avgGood, avgBad, avgSkip, tbl = ScreenComposite.ScreenFromDetails(
          compos, self.details)
        self._dump(1, nGood, misCount, nSkipped, avgGood, avgBad, avgSkip, tbl)
        self.assertEqual(nGood, 92)
        self.assertEqual(misCount, 3)
        self.assertEqual(nSkipped, 0)
        self.assertAlmostEqual(avgGood, 0.9804, 4)
        self.assertAlmostEqual(avgBad, 0.6667, 4)
        self.assertAlmostEqual(avgSkip, 0.0000, 4)
        self.assertEqual(tbl[0, 0], 55)
        self.assertEqual(tbl[1, 1], 37)
        self.assertEqual(tbl[0, 1], 1)
        self.assertEqual(tbl[1, 0], 2)

    def test2_include_holdout(self):
        # """ include holdout data only """
        self.details.tableName = 'ferro_quant'
        self.details.doHoldout = 1
        self.details.doTraining = 0

        compos = self._open(os.path.join(self.baseDir, 'ferromag_quant_10.pkl'))
        tgt = 7
        self.assertEqual(len(compos), tgt)

        nGood, misCount, nSkipped, avgGood, avgBad, avgSkip, tbl = ScreenComposite.ScreenFromDetails(
          compos, self.details)
        self._dump(2, nGood, misCount, nSkipped, avgGood, avgBad, avgSkip, tbl)
        self.assertEqual(nGood, 27)
        self.assertEqual(misCount, 2)
        self.assertEqual(nSkipped, 0)
        self.assertAlmostEqual(avgGood, 0.9778, 4)
        self.assertAlmostEqual(avgBad, 0.7000, 4)
        self.assertAlmostEqual(avgSkip, 0.0000, 4)
        self.assertEqual(tbl[0, 0], 18)
        self.assertEqual(tbl[1, 1], 9)
        self.assertEqual(tbl[0, 1], 1)
        self.assertEqual(tbl[1, 0], 1)

    def test3_include_training(self):
        # """ include training data only """
        self.details.tableName = 'ferro_quant'
        self.details.doHoldout = 0
        self.details.doTraining = 1

        compos = self._open(os.path.join(self.baseDir, 'ferromag_quant_10.pkl'))
        tgt = 7
        self.assertEqual(len(compos), tgt, 'bad composite loaded: %d != %d' % (len(compos), tgt))

        nGood, misCount, nSkipped, avgGood, avgBad, avgSkip, tbl = ScreenComposite.ScreenFromDetails(
          compos, self.details)
        self._dump(3, nGood, misCount, nSkipped, avgGood, avgBad, avgSkip, tbl)
        self.assertEqual(nGood, 65)
        self.assertEqual(misCount, 1)
        self.assertEqual(nSkipped, 0)
        self.assertAlmostEqual(avgGood, 0.9815, 4)
        self.assertAlmostEqual(avgBad, 0.6000, 4)
        self.assertAlmostEqual(avgSkip, 0.0000, 4)
        self.assertEqual(tbl[0, 0], 37)
        self.assertEqual(tbl[1, 1], 28)
        self.assertEqual(tbl[0, 1], 0)
        self.assertEqual(tbl[1, 0], 1)

    def test4_thresholding(self):
        # """ include thresholding """
        self.details.tableName = 'ferro_quant'
        self.details.threshold = 0.80
        self.details.doHoldout = 0
        self.details.doTraining = 0

        compos = self._open(os.path.join(self.baseDir, 'ferromag_quant_10.pkl'))
        tgt = 7
        self.assertEqual(len(compos), tgt)

        nGood, misCount, nSkipped, avgGood, avgBad, avgSkip, tbl = ScreenComposite.ScreenFromDetails(
          compos, self.details)
        self._dump(4, nGood, misCount, nSkipped, avgGood, avgBad, avgSkip, tbl)
        self.assertEqual(nGood, 84)
        self.assertEqual(misCount, 0)
        self.assertEqual(nSkipped, 11)
        self.assertAlmostEqual(avgGood, 1.0000, 4)
        self.assertAlmostEqual(avgBad, 0.0000, 4)
        self.assertAlmostEqual(avgSkip, 0.7455, 4)
        self.assertEqual(tbl[0, 0], 50)
        self.assertEqual(tbl[1, 1], 34)
        self.assertEqual(tbl[0, 1], 0)
        self.assertEqual(tbl[1, 0], 0)

    def test5_basics(self):
        # """ basics """
        self.details.tableName = 'ferro_noquant'

        compos = self._open(os.path.join(self.baseDir, 'ferromag_auto_10_3.pkl'))
        tgt = 10
        self.assertEqual(len(compos), tgt)

        tpl = ScreenComposite.ScreenFromDetails(compos, self.details)
        nGood, misCount, nSkipped, avgGood, avgBad, avgSkip, tbl = tpl
        self._dump(5, nGood, misCount, nSkipped, avgGood, avgBad, avgSkip, tbl)
        self.assertEqual(nGood, 99)
        self.assertEqual(misCount, 4)
        self.assertEqual(nSkipped, 0)
        self.assertAlmostEqual(avgGood, 0.9485, 4)
        self.assertAlmostEqual(avgBad, 0.7500, 4)
        self.assertAlmostEqual(avgSkip, 0.0000, 4)
        self.assertEqual(tbl[0, 0], 53)
        self.assertEqual(tbl[1, 1], 46)
        self.assertEqual(tbl[0, 1], 2)
        self.assertEqual(tbl[1, 0], 2)

    def test6_multiple_models(self):
        # """ multiple models """
        self.details.tableName = 'ferro_noquant'
        compos = self._open(os.path.join(self.baseDir, 'ferromag_auto_10_3.pkl'))
        tgt = 10
        self.assertEqual(len(compos), tgt)
        composites = [compos, compos]
        tpl = ScreenComposite.ScreenFromDetails(composites, self.details)
        nGood, misCount, nSkipped, avgGood, avgBad, avgSkip, tbl = tpl
        self.assertEqual(nGood[0], 99)
        self.assertEqual(misCount[0], 4)
        self.assertEqual(nSkipped[0], 0)
        self.assertAlmostEqual(avgGood[0], .94848, 4)
        self.assertAlmostEqual(avgBad[0], .75, 4)
        self.assertAlmostEqual(avgSkip[0], 0.0, 4)
        self.assertEqual(nGood[1], 0)
        self.assertEqual(misCount[1], 0)
        self.assertEqual(nSkipped[1], 0)
        self.assertEqual(avgGood[1], 0)
        self.assertEqual(avgBad[1], 0)
        self.assertEqual(avgSkip[1], 0)
        self.assertEqual(tbl[0, 0], 53)
        self.assertEqual(tbl[1, 1], 46)
        self.assertEqual(tbl[0, 1], 2)
        self.assertEqual(tbl[1, 0], 2)

    def test7_shuffle(self):
        # """ shuffle """
        self.details.tableName = 'ferro_noquant'
        compos = self._open(os.path.join(self.baseDir, 'ferromag_shuffle_10_3.pkl'))
        tgt = 10
        self.assertEqual(len(compos), tgt)
        self.details.shuffleActivities = 1
        nGood, misCount, nSkipped, avgGood, avgBad, avgSkip, tbl = ScreenComposite.ScreenFromDetails(
          compos, self.details)
        self._dump(7, nGood, misCount, nSkipped, avgGood, avgBad, avgSkip, tbl)
        self.assertEqual(nGood, 52)
        self.assertEqual(misCount, 51)
        self.assertEqual(nSkipped, 0)
        self.assertAlmostEqual(avgGood, 0.7308, 4)
        self.assertAlmostEqual(avgBad, 0.7745, 4)
        self.assertAlmostEqual(avgSkip, 0.0000, 4)
        self.assertEqual(tbl[0, 0], 31)
        self.assertEqual(tbl[1, 1], 21)
        self.assertEqual(tbl[0, 1], 24)
        self.assertEqual(tbl[1, 0], 27)

    def test8_shuffle_segmentation(self):
        # """ shuffle with segmentation """
        self.details.tableName = 'ferro_noquant'
        compos = self._open(os.path.join(self.baseDir, 'ferromag_shuffle_10_3.pkl'))
        tgt = 10
        self.assertEqual(len(compos), tgt)
        self.details.shuffleActivities = 1
        self.details.doHoldout = 1
        nGood, misCount, nSkipped, avgGood, avgBad, avgSkip, tbl = ScreenComposite.ScreenFromDetails(
          compos, self.details)
        self._dump(8, nGood, misCount, nSkipped, avgGood, avgBad, avgSkip, tbl)
        self.assertEqual(nGood, 18)
        self.assertEqual(misCount, 13)
        self.assertEqual(nSkipped, 0)
        self.assertAlmostEqual(avgGood, 0.7444, 4)
        self.assertAlmostEqual(avgBad, 0.8538, 4)
        self.assertAlmostEqual(avgSkip, 0.0000, 4)
        self.assertEqual(tbl[0, 0], 11)
        self.assertEqual(tbl[1, 1], 7)
        self.assertEqual(tbl[0, 1], 5)
        self.assertEqual(tbl[1, 0], 8)

    def test9_shuffle_segmentation2(self):
        # """ shuffle with segmentation2 """
        self.details.tableName = 'ferro_noquant'
        compos = self._open(os.path.join(self.baseDir, 'ferromag_shuffle_10_3.pkl'))
        tgt = 10
        self.assertEqual(len(compos), tgt)
        self.details.shuffleActivities = 1
        self.details.doTraining = 1
        nGood, misCount, nSkipped, avgGood, avgBad, avgSkip, tbl = ScreenComposite.ScreenFromDetails(
          compos, self.details)
        self._dump(9, nGood, misCount, nSkipped, avgGood, avgBad, avgSkip, tbl)
        self.assertEqual(nGood, 34)
        self.assertEqual(misCount, 38)
        self.assertEqual(nSkipped, 0)
        self.assertAlmostEqual(avgGood, 0.7235, 4)
        self.assertAlmostEqual(avgBad, 0.7474, 4)
        self.assertAlmostEqual(avgSkip, 0.0000, 4)
        self.assertEqual(tbl[0, 0], 20)
        self.assertEqual(tbl[1, 1], 14)
        self.assertEqual(tbl[0, 1], 19)
        self.assertEqual(tbl[1, 0], 19)

    def test10_filtering(self):
        # """ filtering """
        self.details.tableName = 'ferro_noquant'
        compos = self._open(os.path.join(self.baseDir, 'ferromag_filt_10_3.pkl'))
        tgt = 10
        self.assertEqual(len(compos), tgt)
        self.details.filterVal = 1
        self.details.filterFrac = .33

        nGood, misCount, nSkipped, avgGood, avgBad, avgSkip, tbl = ScreenComposite.ScreenFromDetails(
          compos, self.details)
        self._dump(10, nGood, misCount, nSkipped, avgGood, avgBad, avgSkip, tbl)
        self.assertEqual(nGood, 90)
        self.assertEqual(misCount, 13)
        self.assertEqual(nSkipped, 0)
        self.assertAlmostEqual(avgGood, 0.9578, 4)
        self.assertAlmostEqual(avgBad, 0.8538, 4)
        self.assertAlmostEqual(avgSkip, 0.0000, 4)
        self.assertEqual(tbl[0, 0], 54)
        self.assertEqual(tbl[1, 1], 36)
        self.assertEqual(tbl[0, 1], 1)
        self.assertEqual(tbl[1, 0], 12)

    def test11_filtering_segmentation(self):
        # """ filtering with segmentation """
        self.details.tableName = 'ferro_noquant'
        compos = self._open(os.path.join(self.baseDir, 'ferromag_filt_10_3.pkl'))
        tgt = 10
        self.assertEqual(len(compos), tgt)
        self.details.doHoldout = 1
        self.details.filterVal = 1
        self.details.filterFrac = .33

        nGood, misCount, nSkipped, avgGood, avgBad, avgSkip, tbl = ScreenComposite.ScreenFromDetails(
          compos, self.details)
        self._dump(11, nGood, misCount, nSkipped, avgGood, avgBad, avgSkip, tbl)
        self.assertEqual(nGood, 40)
        self.assertEqual(misCount, 6)
        self.assertEqual(nSkipped, 0)
        self.assertAlmostEqual(avgGood, 0.9550, 4)
        self.assertAlmostEqual(avgBad, 0.8833, 4)
        self.assertAlmostEqual(avgSkip, 0.0000, 4)
        self.assertEqual(tbl[0, 0], 17)
        self.assertEqual(tbl[1, 1], 23)
        self.assertEqual(tbl[0, 1], 0)
        self.assertEqual(tbl[1, 0], 6)

    def test12_naiveBayes_composite(self):
        # """ test the naive bayes composite"""
        self.details.tableName = 'ferro_noquant'
        compos = self._open(os.path.join(self.baseDir, 'ferromag_NaiveBayes.pkl'))
        tgt = 10
        self.assertEqual(len(compos), tgt)
        self.details.doHoldout = 1
        nGood, misCount, nSkipped, avgGood, avgBad, avgSkip, tbl = ScreenComposite.ScreenFromDetails(
          compos, self.details)
        self._dump(12, nGood, misCount, nSkipped, avgGood, avgBad, avgSkip, tbl)
        self.assertEqual(nGood, 25)
        self.assertEqual(misCount, 6)
        self.assertEqual(nSkipped, 0)
        self.assertAlmostEqual(avgGood, 0.9560, 4)
        self.assertAlmostEqual(avgBad, 0.8333, 4)
        self.assertAlmostEqual(avgSkip, 0.0000, 4)
        self.assertEqual(tbl[0, 0], 14)
        self.assertEqual(tbl[1, 1], 11)
        self.assertEqual(tbl[0, 1], 6)
        self.assertEqual(tbl[1, 0], 0)


if __name__ == '__main__':  # pragma: nocover
    unittest.main()
