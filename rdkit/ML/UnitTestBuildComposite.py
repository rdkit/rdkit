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
"""unit testing code for the BuildComposite functionality

"""
import io
import os
import unittest

from rdkit import RDConfig
from rdkit.Dbase.DbConnection import DbConnect
from rdkit.ML import BuildComposite
import pickle


class TestCase(unittest.TestCase):

    def setUp(self):
        self.baseDir = os.path.join(RDConfig.RDCodeDir, 'ML', 'test_data')
        self.dbName = RDConfig.RDTestDatabase

        self.details = BuildComposite.SetDefaults()
        self.details.dbName = self.dbName
        self.details.dbUser = RDConfig.defaultDBUser
        self.details.dbPassword = RDConfig.defaultDBPassword

    def _init(self, refCompos, copyBounds=0):
        BuildComposite._verbose = 0
        conn = DbConnect(self.details.dbName, self.details.tableName)
        cols = [x.upper() for x in conn.GetColumnNames()]
        cDescs = [x.upper() for x in refCompos.GetDescriptorNames()]
        self.assertEqual(cols, cDescs)

        self.details.nModels = 10
        self.details.lockRandom = 1
        self.details.randomSeed = refCompos._randomSeed
        self.details.splitFrac = refCompos._splitFrac

        if self.details.splitFrac:
            self.details.splitRun = 1
        else:
            self.details.splitRun = 0

        if not copyBounds:
            self.details.qBounds = [0] * len(cols)
        else:
            self.details.qBounds = refCompos.GetQuantBounds()[0]

    def compare(self, compos, refCompos):
        self.assertEqual(len(compos), len(refCompos))
        cs = []
        rcs = []
        for i in range(len(compos)):
            cs.append(compos[i])
            rcs.append(refCompos[i])

        cs.sort(key=lambda x: (x[2], x[2]))
        rcs.sort(key=lambda x: (x[2], x[2]))

        for i in range(len(compos)):
            _, count, err = cs[i]
            _, refCount, refErr = rcs[i]
            self.assertEqual(count, refCount)
            self.assertAlmostEqual(err, refErr, 4)

    def test1_basics(self):
        # """ basics """
        self.details.tableName = 'ferro_quant'
        refComposName = 'ferromag_quant_10.pkl'

        with open(os.path.join(self.baseDir, refComposName), 'r') as pklTF:
            buf = pklTF.read().replace('\r\n', '\n').encode('utf-8')
            pklTF.close()
        with io.BytesIO(buf) as pklF:
            refCompos = pickle.load(pklF)

        # first make sure the data are intact
        self._init(refCompos)
        compos = BuildComposite.RunIt(self.details, saveIt=0)

        # pickle.dump(compos,open(os.path.join(self.baseDir,refComposName), 'wb'))
        # with open(os.path.join(self.baseDir,refComposName), 'rb') as pklF:
        #   refCompos = pickle.load(pklF)

        self.compare(compos, refCompos)

    def test2_depth_limit(self):
        # """ depth limit """
        self.details.tableName = 'ferro_quant'
        refComposName = 'ferromag_quant_10_3.pkl'

        with open(os.path.join(self.baseDir, refComposName), 'r') as pklTF:
            buf = pklTF.read().replace('\r\n', '\n').encode('utf-8')
            pklTF.close()
        with io.BytesIO(buf) as pklF:
            refCompos = pickle.load(pklF)

        # first make sure the data are intact
        self._init(refCompos)
        self.details.limitDepth = 3
        compos = BuildComposite.RunIt(self.details, saveIt=0)

        self.compare(compos, refCompos)

    def test3_depth_limit_less_greedy(self):
        # """ depth limit + less greedy """
        self.details.tableName = 'ferro_quant'
        refComposName = 'ferromag_quant_10_3_lessgreedy.pkl'

        with open(os.path.join(self.baseDir, refComposName), 'r') as pklTF:
            buf = pklTF.read().replace('\r\n', '\n').encode('utf-8')
            pklTF.close()
        with io.BytesIO(buf) as pklF:
            refCompos = pickle.load(pklF)

        # first make sure the data are intact
        self._init(refCompos)
        self.details.limitDepth = 3
        self.details.lessGreedy = 1
        compos = BuildComposite.RunIt(self.details, saveIt=0)

        self.compare(compos, refCompos)

    def test4_more_trees(self):
        # """ more trees """
        self.details.tableName = 'ferro_quant'
        refComposName = 'ferromag_quant_50_3.pkl'

        with open(os.path.join(self.baseDir, refComposName), 'r') as pklTF:
            buf = pklTF.read().replace('\r\n', '\n').encode('utf-8')
            pklTF.close()
        with io.BytesIO(buf) as pklF:
            refCompos = pickle.load(pklF)

        # first make sure the data are intact
        self._init(refCompos)
        self.details.limitDepth = 3
        self.details.nModels = 50
        compos = BuildComposite.RunIt(self.details, saveIt=0)

        self.compare(compos, refCompos)

    def test5_auto_bounds(self):
        # """ auto bounds """
        self.details.tableName = 'ferro_noquant'
        refComposName = 'ferromag_auto_10_3.pkl'

        with open(os.path.join(self.baseDir, refComposName), 'r') as pklTF:
            buf = pklTF.read().replace('\r\n', '\n').encode('utf-8')
            pklTF.close()
        with io.BytesIO(buf) as pklF:
            refCompos = pickle.load(pklF)

        # first make sure the data are intact
        self._init(refCompos, copyBounds=1)
        self.details.limitDepth = 3
        self.details.nModels = 10
        compos = BuildComposite.RunIt(self.details, saveIt=0)

        self.compare(compos, refCompos)

    def test6_auto_bounds_real_activity(self):
        # """ auto bounds with a real valued activity"""
        self.details.tableName = 'ferro_noquant_realact'
        refComposName = 'ferromag_auto_10_3.pkl'

        with open(os.path.join(self.baseDir, refComposName), 'r') as pklTF:
            buf = pklTF.read().replace('\r\n', '\n').encode('utf-8')
            pklTF.close()
        with io.BytesIO(buf) as pklF:
            refCompos = pickle.load(pklF)

        # first make sure the data are intact
        self._init(refCompos, copyBounds=1)
        self.details.limitDepth = 3
        self.details.nModels = 10
        self.details.activityBounds = [0.5]
        compos = BuildComposite.RunIt(self.details, saveIt=0)

        self.compare(compos, refCompos)

    def test7_composite_naiveBayes(self):
        # """ Test composite of naive bayes"""
        self.details.tableName = 'ferro_noquant'
        refComposName = 'ferromag_NaiveBayes.pkl'
        with open(os.path.join(self.baseDir, refComposName), 'r') as pklTFile:
            buf = pklTFile.read().replace('\r\n', '\n').encode('utf-8')
            pklTFile.close()
        with io.BytesIO(buf) as pklFile:
            refCompos = pickle.load(pklFile)
        self._init(refCompos, copyBounds=1)
        self.details.useTrees = 0
        self.details.useNaiveBayes = 1
        self.details.mEstimateVal = 20.0
        self.details.qBounds = [0] + [2] * 6 + [0]
        compos = BuildComposite.RunIt(self.details, saveIt=0)

        self.compare(compos, refCompos)


if __name__ == '__main__':  # pragma: nocover
    unittest.main()
