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


class TestCase(unittest.TestCase):

    def setUp(self):
        self.baseDir = os.path.join(RDConfig.RDCodeDir, 'ML', 'test_data')
        self.dbName = RDConfig.RDTestDatabase
        self.details = ScreenComposite.SetDefaults()
        self.details.dbName = self.dbName
        self.details.dbUser = RDConfig.defaultDBUser
        self.details.dbPassword = RDConfig.defaultDBPassword

    def test1_basics(self):
        # """ basics """
        self.details.tableName = 'ferro_quant'
        with open(os.path.join(self.baseDir, 'ferromag_quant_10.pkl'), 'r') as pklTF:
            buf = pklTF.read().replace('\r\n', '\n').encode('utf-8')
            pklTF.close()
        with io.BytesIO(buf) as pklF:
            compos = pickle.load(pklF)
        tgt = 5
        self.assertEqual(len(compos), tgt)

        nGood, misCount, nSkipped, avgGood, avgBad, avgSkip, tbl = ScreenComposite.ScreenFromDetails(
          compos, self.details)
        self.assertEqual(nGood, 93)
        self.assertEqual(misCount, 2)
        self.assertEqual(nSkipped, 0)
        self.assertAlmostEqual(avgGood, .9871, 4)
        self.assertAlmostEqual(avgBad, .8000, 4)
        self.assertAlmostEqual(avgSkip, 0, 4)
        self.assertEqual(tbl[0, 0], 54)
        self.assertEqual(tbl[1, 1], 39)
        self.assertEqual(tbl[0, 1], 2)
        self.assertEqual(tbl[1, 0], 0)

    def test2_include_holdout(self):
        # """ include holdout data only """
        self.details.tableName = 'ferro_quant'
        self.details.doHoldout = 1
        self.details.doTraining = 0

        with open(os.path.join(self.baseDir, 'ferromag_quant_10.pkl'), 'r') as pklTF:
            buf = pklTF.read().replace('\r\n', '\n').encode('utf-8')
            pklTF.close()
        with io.BytesIO(buf) as pklF:
            compos = pickle.load(pklF)
        tgt = 5
        self.assertEqual(len(compos), tgt)

        nGood, misCount, nSkipped, avgGood, avgBad, avgSkip, tbl = ScreenComposite.ScreenFromDetails(
          compos, self.details)
        self.assertEqual(nGood, 28)
        self.assertEqual(misCount, 1)
        self.assertEqual(nSkipped, 0)
        self.assertAlmostEqual(avgGood, .9964, 4)
        self.assertAlmostEqual(avgBad, 1.000, 4)
        self.assertAlmostEqual(avgSkip, 0, 4)
        self.assertEqual(tbl[0, 0], 16)
        self.assertEqual(tbl[1, 1], 12)
        self.assertEqual(tbl[0, 1], 1)
        self.assertEqual(tbl[1, 0], 0)

    def test3_include_training(self):
        # """ include training data only """
        self.details.tableName = 'ferro_quant'
        self.details.doHoldout = 0
        self.details.doTraining = 1

        with open(os.path.join(self.baseDir, 'ferromag_quant_10.pkl'), 'r') as pklTF:
            buf = pklTF.read().replace('\r\n', '\n').encode('utf-8')
            pklTF.close()
        with io.BytesIO(buf) as pklF:
            compos = pickle.load(pklF)
        tgt = 5
        self.assertEqual(len(compos), tgt, 'bad composite loaded: %d != %d' % (len(compos), tgt))

        nGood, misCount, nSkipped, avgGood, avgBad, avgSkip, tbl = ScreenComposite.ScreenFromDetails(
          compos, self.details)
        self.assertEqual(nGood, 65)
        self.assertEqual(misCount, 1)
        self.assertEqual(nSkipped, 0)
        self.assertAlmostEqual(avgGood, .98307, 4)
        self.assertAlmostEqual(avgBad, 0.600, 4)
        self.assertAlmostEqual(avgSkip, 0, 4)
        self.assertEqual(tbl[0, 0], 38, tbl)
        self.assertEqual(tbl[1, 1], 27)
        self.assertEqual(tbl[0, 1], 1)
        self.assertEqual(tbl[1, 0], 0)

    def test4_thresholding(self):
        # """ include thresholding """
        self.details.tableName = 'ferro_quant'
        self.details.threshold = 0.80
        self.details.doHoldout = 0
        self.details.doTraining = 0

        with open(os.path.join(self.baseDir, 'ferromag_quant_10.pkl'), 'r') as pklTF:
            buf = pklTF.read().replace('\r\n', '\n').encode('utf-8')
            pklTF.close()
        with io.BytesIO(buf) as pklF:
            compos = pickle.load(pklF)
        tgt = 5
        self.assertEqual(len(compos), tgt)

        nGood, misCount, nSkipped, avgGood, avgBad, avgSkip, tbl = ScreenComposite.ScreenFromDetails(
          compos, self.details)
        self.assertEqual(nGood, 91)
        self.assertEqual(misCount, 1)
        self.assertEqual(nSkipped, 3)
        self.assertAlmostEqual(avgGood, 0.9956, 4)
        self.assertAlmostEqual(avgBad, 1.000, 4)
        self.assertAlmostEqual(avgSkip, 0.6000, 4)
        self.assertEqual(tbl[0, 0], 54)
        self.assertEqual(tbl[1, 1], 37)
        self.assertEqual(tbl[0, 1], 1)
        self.assertEqual(tbl[1, 0], 0)

    def test5_basics(self):
        # """ basics """
        self.details.tableName = 'ferro_noquant'

        with open(os.path.join(self.baseDir, 'ferromag_auto_10_3.pkl'), 'r') as pklTF:
            buf = pklTF.read().replace('\r\n', '\n').encode('utf-8')
            pklTF.close()
        with io.BytesIO(buf) as pklF:
            compos = pickle.load(pklF)
        tgt = 10
        self.assertEqual(len(compos), tgt)

        tpl = ScreenComposite.ScreenFromDetails(compos, self.details)
        nGood, misCount, nSkipped, avgGood, avgBad, avgSkip, tbl = tpl

        self.assertEqual(nGood, 95)
        self.assertEqual(misCount, 8)
        self.assertEqual(nSkipped, 0)
        self.assertAlmostEqual(avgGood, .9684, 4)
        self.assertAlmostEqual(avgBad, .8375, 4)
        self.assertAlmostEqual(avgSkip, 0, 4)
        self.assertEqual(tbl[0, 0], 50)
        self.assertEqual(tbl[1, 1], 45)
        self.assertEqual(tbl[0, 1], 5)
        self.assertEqual(tbl[1, 0], 3)

    def test6_multiple_models(self):
        # """ multiple models """
        self.details.tableName = 'ferro_noquant'
        with open(os.path.join(self.baseDir, 'ferromag_auto_10_3.pkl'), 'r') as pklTF:
            buf = pklTF.read().replace('\r\n', '\n').encode('utf-8')
            pklTF.close()
        with io.BytesIO(buf) as pklF:
            compos = pickle.load(pklF)
        tgt = 10
        self.assertEqual(len(compos), tgt)
        composites = [compos, compos]
        tpl = ScreenComposite.ScreenFromDetails(composites, self.details)
        nGood, misCount, nSkipped, avgGood, avgBad, avgSkip, tbl = tpl
        self.assertEqual(nGood[0], 95)
        self.assertEqual(misCount[0], 8)
        self.assertEqual(nSkipped[0], 0)
        self.assertAlmostEqual(avgGood[0], .9684, 4)
        self.assertAlmostEqual(avgBad[0], .8375, 4)
        self.assertAlmostEqual(avgSkip[0], 0.0, 4)
        self.assertEqual(nGood[1], 0)
        self.assertEqual(misCount[1], 0)
        self.assertEqual(nSkipped[1], 0)
        self.assertEqual(avgGood[1], 0)
        self.assertEqual(avgBad[1], 0)
        self.assertEqual(avgSkip[1], 0)
        self.assertEqual(tbl[0, 0], 50)
        self.assertEqual(tbl[1, 1], 45)
        self.assertEqual(tbl[0, 1], 5)
        self.assertEqual(tbl[1, 0], 3)

    def test7_shuffle(self):
        # """ shuffle """
        self.details.tableName = 'ferro_noquant'
        with open(os.path.join(self.baseDir, 'ferromag_shuffle_10_3.pkl'), 'r') as pklTF:
            buf = pklTF.read().replace('\r\n', '\n').encode('utf-8')
            pklTF.close()
        with io.BytesIO(buf) as pklF:
            compos = pickle.load(pklF)
        tgt = 10
        self.assertEqual(len(compos), tgt)
        self.details.shuffleActivities = 1
        nGood, misCount, nSkipped, avgGood, avgBad, avgSkip, tbl = ScreenComposite.ScreenFromDetails(
          compos, self.details)
        self.assertEqual(nGood, 50)
        self.assertEqual(misCount, 53)
        self.assertEqual(nSkipped, 0)
        self.assertAlmostEqual(avgGood, .7380, 4)
        self.assertAlmostEqual(avgBad, .7660, 4)
        self.assertAlmostEqual(avgSkip, 0, 4)
        self.assertEqual(tbl[0, 0], 30)
        self.assertEqual(tbl[1, 1], 20)
        self.assertEqual(tbl[0, 1], 25)
        self.assertEqual(tbl[1, 0], 28)

    def test8_shuffle_segmentation(self):
        # """ shuffle with segmentation """
        self.details.tableName = 'ferro_noquant'
        with open(os.path.join(self.baseDir, 'ferromag_shuffle_10_3.pkl'), 'r') as pklTF:
            buf = pklTF.read().replace('\r\n', '\n').encode('utf-8')
            pklTF.close()
        with io.BytesIO(buf) as pklF:
            compos = pickle.load(pklF)
        tgt = 10
        self.assertEqual(len(compos), tgt)
        self.details.shuffleActivities = 1
        self.details.doHoldout = 1
        nGood, misCount, nSkipped, avgGood, avgBad, avgSkip, tbl = ScreenComposite.ScreenFromDetails(
          compos, self.details)
        self.assertEqual(nGood, 19)
        self.assertEqual(misCount, 12)
        self.assertEqual(nSkipped, 0)
        self.assertAlmostEqual(avgGood, .7737, 4)
        self.assertAlmostEqual(avgBad, .7500, 4)
        self.assertAlmostEqual(avgSkip, 0, 4)
        self.assertEqual(tbl[0, 0], 12)
        self.assertEqual(tbl[1, 1], 7)
        self.assertEqual(tbl[0, 1], 6)
        self.assertEqual(tbl[1, 0], 6)

    def test9_shuffle_segmentation2(self):
        # """ shuffle with segmentation2 """
        self.details.tableName = 'ferro_noquant'
        with open(os.path.join(self.baseDir, 'ferromag_shuffle_10_3.pkl'), 'r') as pklTF:
            buf = pklTF.read().replace('\r\n', '\n').encode('utf-8')
            pklTF.close()
        with io.BytesIO(buf) as pklF:
            compos = pickle.load(pklF)
        tgt = 10
        self.assertEqual(len(compos), tgt)
        self.details.shuffleActivities = 1
        self.details.doTraining = 1
        nGood, misCount, nSkipped, avgGood, avgBad, avgSkip, tbl = ScreenComposite.ScreenFromDetails(
          compos, self.details)
        self.assertEqual(nGood, 31)
        self.assertEqual(misCount, 41)
        self.assertEqual(nSkipped, 0)
        self.assertAlmostEqual(avgGood, .7161, 4)
        self.assertAlmostEqual(avgBad, .7707, 4)
        self.assertAlmostEqual(avgSkip, 0.0, 4)
        self.assertEqual(tbl[0, 0], 18)
        self.assertEqual(tbl[1, 1], 13)
        self.assertEqual(tbl[0, 1], 19)
        self.assertEqual(tbl[1, 0], 22)

    def test10_filtering(self):
        # """ filtering """
        self.details.tableName = 'ferro_noquant'
        with open(os.path.join(self.baseDir, 'ferromag_filt_10_3.pkl'), 'r') as pklTF:
            buf = pklTF.read().replace('\r\n', '\n').encode('utf-8')
            pklTF.close()
        with io.BytesIO(buf) as pklF:
            compos = pickle.load(pklF)
        tgt = 10
        self.assertEqual(len(compos), tgt)
        self.details.filterVal = 1
        self.details.filterFrac = .33

        nGood, misCount, nSkipped, avgGood, avgBad, avgSkip, tbl = ScreenComposite.ScreenFromDetails(
          compos, self.details)
        self.assertEqual(nGood, 90)
        self.assertEqual(misCount, 13)
        self.assertEqual(nSkipped, 0)
        self.assertAlmostEqual(avgGood, .9578, 4)
        self.assertAlmostEqual(avgBad, .8538, 4)
        self.assertAlmostEqual(avgSkip, 0, 4)
        self.assertEqual(tbl[0, 0], 54)
        self.assertEqual(tbl[1, 1], 36)
        self.assertEqual(tbl[0, 1], 1)
        self.assertEqual(tbl[1, 0], 12)

    def test11_filtering_segmentation(self):
        # """ filtering with segmentation """
        self.details.tableName = 'ferro_noquant'
        with open(os.path.join(self.baseDir, 'ferromag_filt_10_3.pkl'), 'r') as pklTF:
            buf = pklTF.read().replace('\r\n', '\n').encode('utf-8')
            pklTF.close()
        with io.BytesIO(buf) as pklF:
            compos = pickle.load(pklF)
        tgt = 10
        self.assertEqual(len(compos), tgt)
        self.details.doHoldout = 1
        self.details.filterVal = 1
        self.details.filterFrac = .33

        nGood, misCount, nSkipped, avgGood, avgBad, avgSkip, tbl = ScreenComposite.ScreenFromDetails(
          compos, self.details)

        self.assertEqual(nGood, 37)
        self.assertEqual(misCount, 6)
        self.assertEqual(nSkipped, 0)
        self.assertAlmostEqual(avgGood, .95946, 4)
        self.assertAlmostEqual(avgBad, .85, 4)
        self.assertAlmostEqual(avgSkip, 0, 4)
        self.assertEqual(tbl[0, 0], 14)
        self.assertEqual(tbl[1, 1], 23)
        self.assertEqual(tbl[0, 1], 1)
        self.assertEqual(tbl[1, 0], 5)

    def test12_naiveBayes_composite(self):
        # """ test the naive bayes composite"""
        self.details.tableName = 'ferro_noquant'
        with open(os.path.join(self.baseDir, 'ferromag_NaiveBayes.pkl'), 'r') as pklTF:
            buf = pklTF.read().replace('\r\n', '\n').encode('utf-8')
            pklTF.close()
        with io.BytesIO(buf) as pklF:
            compos = pickle.load(pklF)
        tgt = 10
        self.assertEqual(len(compos), tgt)
        self.details.doHoldout = 1
        nGood, misCount, nSkipped, avgGood, avgBad, avgSkip, tbl = ScreenComposite.ScreenFromDetails(
          compos, self.details)
        self.assertEqual(nGood, 25)
        self.assertEqual(misCount, 6)
        self.assertEqual(nSkipped, 0)
        self.assertAlmostEqual(avgGood, 0.9800, 4)
        self.assertAlmostEqual(avgBad, 0.86667, 4)
        self.assertAlmostEqual(avgSkip, 0, 4)
        self.assertEqual(tbl[0, 0], 9)
        self.assertEqual(tbl[0, 1], 6)
        self.assertEqual(tbl[1, 0], 0)
        self.assertEqual(tbl[1, 1], 16)


if __name__ == '__main__':  # pragma: nocover
    unittest.main()
