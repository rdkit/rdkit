# $Id$
#
# Copyright (C) 2003-2006 greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" unit testing code for BitEnsembles
"""
import os
import shutil
import tempfile
import unittest

from rdkit import RDConfig
from rdkit.DataStructs import SparseBitVect
# This import is important to initialize the BitEnsemble module
from rdkit.DataStructs import BitEnsembleDb
from rdkit.DataStructs.BitEnsemble import BitEnsemble


class TestCase(unittest.TestCase):

  def test1(self):
    ensemble = BitEnsemble()
    ensemble.SetBits([1, 11, 21, 31])
    self.assertEqual(ensemble.GetNumBits(), 4)
    bv = SparseBitVect(100)
    bv.SetBit(1)
    bv.SetBit(11)
    bv.SetBit(13)

    score = ensemble.ScoreWithOnBits(bv)
    assert score == 2, 'bad score: %d' % (score)
    score = ensemble.ScoreWithIndex(bv)
    assert score == 2, 'bad score: %d' % (score)

  def test2(self):
    ensemble = BitEnsemble([1, 11, 21, 31])
    bv = SparseBitVect(100)
    bv.SetBit(1)
    bv.SetBit(11)
    bv.SetBit(13)

    score = ensemble.ScoreWithOnBits(bv)
    assert score == 2, 'bad score: %d' % (score)
    score = ensemble.ScoreWithIndex(bv)
    assert score == 2, 'bad score: %d' % (score)

  def test3(self):
    ensemble = BitEnsemble()
    for bit in [1, 11, 21, 31]:
      ensemble.AddBit(bit)
    bv = SparseBitVect(100)
    bv.SetBit(1)
    bv.SetBit(11)
    bv.SetBit(13)

    score = ensemble.ScoreWithOnBits(bv)
    assert score == 2, 'bad score: %d' % (score)
    score = ensemble.ScoreWithIndex(bv)
    assert score == 2, 'bad score: %d' % (score)

  def _setupDb(self):
    from rdkit.Dbase.DbConnection import DbConnect
    fName = RDConfig.RDTestDatabase
    if RDConfig.useSqlLite:
      _, tempName = tempfile.mkstemp(suffix='sqlt')
      self.tempDbName = tempName
      shutil.copyfile(fName, tempName)
    else:  # pragma: nocover
      tempName = '::RDTests'
    self.conn = DbConnect(tempName)
    self.dbTblName = 'bit_ensemble_test'
    return self.conn

  def tearDown(self):
    if hasattr(self, 'tempDbName') and RDConfig.useSqlLite and os.path.exists(self.tempDbName):
      try:
        os.unlink(self.tempDbName)
      except Exception:  # pragma: nocover
        import traceback
        traceback.print_exc()

  def testdb1(self):
    """ test the sig - db functionality """
    conn = self._setupDb()
    ensemble = BitEnsemble()
    for bit in [1, 3, 4]:
      ensemble.AddBit(bit)

    sigBs = [([0, 0, 0, 0, 0, 0], (0, 0, 0)),
             ([0, 1, 0, 1, 0, 0], (1, 1, 0)),
             ([0, 1, 0, 0, 1, 0], (1, 0, 1)),
             ([0, 1, 0, 0, 1, 1], (1, 0, 1)), ]
    ensemble.InitScoreTable(conn, self.dbTblName)
    for bs, tgt in sigBs:
      ensemble.ScoreToDb(bs, conn)
    conn.Commit()

    d = conn.GetData(table=self.dbTblName)
    assert len(d) == len(sigBs), 'bad number of results returned'
    for i in range(len(sigBs)):
      bs, tgt = tuple(sigBs[i])
      dbRes = tuple(d[i])
      assert dbRes == tgt, 'bad bits returned: %s != %s' % (str(dbRes), str(tgt))
    d = None
    self.conn = None

  def testdb2(self):
    """ test the sig - db functionality """
    conn = self._setupDb()
    ensemble = BitEnsemble()
    for bit in [1, 3, 4]:
      ensemble.AddBit(bit)

    sigBs = [([0, 0, 0, 0, 0, 0], (0, 0, 0)),
             ([0, 1, 0, 1, 0, 0], (1, 1, 0)),
             ([0, 1, 0, 0, 1, 0], (1, 0, 1)),
             ([0, 1, 0, 0, 1, 1], (1, 0, 1)), ]
    ensemble.InitScoreTable(conn, self.dbTblName, idInfo='id varchar(10)', actInfo='act int')
    for bs, tgt in sigBs:
      ensemble.ScoreToDb(bs, conn, id='foo', act=1)
    conn.Commit()

    d = conn.GetData(table=self.dbTblName)
    assert len(d) == len(sigBs), 'bad number of results returned'
    for i in range(len(sigBs)):
      bs, tgt = tuple(sigBs[i])
      dbRes = tuple(d[i])
      assert dbRes[1:-1] == tgt, 'bad bits returned: %s != %s' % (str(dbRes[1:-1]), str(tgt))

    d = None
    self.conn = None


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
