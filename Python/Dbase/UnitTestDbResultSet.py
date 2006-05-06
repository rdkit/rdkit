# $Id: UnitTestDbResultSet.py 5010 2006-02-22 15:23:39Z glandrum $
#
#  Copyright (C) 2003-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
"""unit testing code for the DbResultSet object

"""
import RDConfig
import unittest,os
from Dbase.DbConnection import DbConnect
from Dbase.DbResultSet import DbResultSet,RandomAccessDbResultSet

class TestCase(unittest.TestCase):
  def setUp(self):
    if not RDConfig.usePgSQL:
      self.dbName = os.path.join(RDConfig.RDCodeDir,'Dbase','testData','TEST.GDB')
    else:
      self.dbName = "::RDTests"
    self.conn = DbConnect(self.dbName)
    self.curs = self.conn.GetCursor()

  def test1(self):
    """ test indexing in, ensure acceptable error conditions
    """
    cmd = 'select * from ten_elements'
    set = RandomAccessDbResultSet(self.curs,self.conn,cmd)
    for i in range(12):
      try:
        val = set[i]
      except IndexError:
        assert i >= 10

  def test2(self):
    """ 
    """
    cmd = 'select * from ten_elements'
    set = RandomAccessDbResultSet(self.curs,self.conn,cmd)
    assert len(set)==10
    for i in range(len(set)):
      try:
        val = set[i]
      except:
        assert 0

  def test3(self):
    """ 
    """
    cmd = 'select * from ten_elements'
    set = DbResultSet(self.curs,self.conn,cmd)
    r = []
    for thing in set:
      r.append(thing)
    assert len(r)==10
        
  def test4(self):
    """ 
    """
    cmd = 'select * from ten_elements_dups'
    set = DbResultSet(self.curs,self.conn,cmd,removeDups=0)
    r = []
    for thing in set:
      r.append(thing)
    assert len(r)==10
        
  def test5(self):
    """ 
    """
    cmd='select * from ten_elements_dups'
    set = RandomAccessDbResultSet(self.curs,self.conn,cmd,removeDups=0)
    assert len(set)==10
    for i in range(len(set)):
      try:
        val = set[i]
      except:
        assert 0

  def test6(self):
    """ 
    """
    cmd = 'select * from ten_elements_dups'
    set = DbResultSet(self.curs,self.conn,cmd,removeDups=0)
    r = []
    for thing in set:
      r.append(thing)
    assert len(r)==10


if __name__ == '__main__':
  unittest.main()

