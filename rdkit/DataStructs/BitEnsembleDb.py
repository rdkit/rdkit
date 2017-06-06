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
""" This functionality gets mixed into the BitEnsemble class

"""
from rdkit.DataStructs.BitEnsemble import BitEnsemble


def _InitScoreTable(self, dbConn, tableName, idInfo='', actInfo=''):
  """ inializes a db table to store our scores

    idInfo and actInfo should be strings with the definitions of the id and
    activity columns of the table (when desired)

  """
  if idInfo:
    cols = [idInfo]
  else:
    cols = []
  for bit in self.GetBits():
    cols.append('Bit_%d smallint' % (bit))
  if actInfo:
    cols.append(actInfo)
  dbConn.AddTable(tableName, ','.join(cols))
  self._dbTableName = tableName


def _ScoreToDb(self, sig, dbConn, tableName=None, id=None, act=None):
  """ scores the "signature" that is passed in and puts the
  results in the db table

  """
  if tableName is None:
    try:
      tableName = self._dbTableName
    except AttributeError:
      raise ValueError('table name not set in BitEnsemble pre call to ScoreToDb()')
  if id is not None:
    cols = [id]
  else:
    cols = []
  score = 0
  for bit in self.GetBits():
    b = sig[bit]
    cols.append(b)
    score += b
  if act is not None:
    cols.append(act)
  dbConn.InsertData(tableName, cols)


BitEnsemble.InitScoreTable = _InitScoreTable
BitEnsemble.ScoreToDb = _ScoreToDb
