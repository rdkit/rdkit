# $Id$
#
#  Copyright (C) 2003-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""unit testing code for database connection objects

"""
from rdkit import RDConfig
import unittest,os,sys
from rdkit.Dbase.DbConnection import DbConnect

class TestCase(unittest.TestCase):
  def setUp(self):
    #print '\n%s: '%self.shortDescription(),
    #sys.stderr.write('\n%s: \n'%self.shortDescription())
    #sys.stderr.flush()
    self.baseDir = os.path.join(RDConfig.RDCodeDir,'Dbase','test_data')
    self.dbName = RDConfig.RDTestDatabase
    self.colHeads=('int_col','floatCol','strCol')
    self.colTypes=('integer','float','string')

  def testAddTable(self):
    """ tests AddTable and GetTableNames functionalities """
    newTblName = 'NEW_TABLE'
    conn = DbConnect(self.dbName)
    try:
      conn.GetCursor().execute('drop table %s'%(newTblName))
    except:
      pass
    conn.Commit()
    conn.AddTable(newTblName,'id int')
    names = [x.strip() for x in conn.GetTableNames()]
    assert newTblName in names,'name (%s) not found in %s'%(newTblName,str(names))
    conn.GetCursor().execute('drop table %s'%(newTblName))

  def testCursor(self):
    """ tests GetCursor and GetTableNames functionalities """

    viewName = 'TEST_VIEW'
    conn = DbConnect(self.dbName)
    curs = conn.GetCursor()
    assert curs
    try:
      curs.execute('drop view %s'%(viewName))
    except:
      pass
    try:
      curs.execute('create view %s as select val,id from ten_elements'%(viewName))
    except:
      import traceback
      traceback.print_exc()
      assert 0
    conn.Commit()
      
    names = [x.strip() for x in conn.GetTableNames(includeViews=0)]
    assert viewName not in names,'improper view found'
    names = [x.strip() for x in conn.GetTableNames(includeViews=1)]
    assert viewName in names,'improper view found in %s'%(str(names))
    try:
      curs.execute('drop view %s'%(viewName))
    except:
      assert 0,'drop table failed'

  def testGetData1(self):
    """ basic functionality
    """
    conn = DbConnect(self.dbName,'ten_elements')
    d = conn.GetData(randomAccess=1)
    assert len(d)==10
    assert tuple(d[0])==(0,11)
    assert tuple(d[2])==(4,31)
    self.failUnlessRaises(IndexError,lambda:d[11])

  def testGetData2(self):
    """ using removeDups
    """
    conn = DbConnect(self.dbName,'ten_elements_dups')
    d = conn.GetData(randomAccess=1,removeDups=1)
    assert tuple(d[0])==(0,11)
    assert tuple(d[2])==(4,31)
    assert len(d)==10
    self.failUnlessRaises(IndexError,lambda:d[11])
    
  def testGetData3(self):
    """ without removeDups
    """
    conn = DbConnect(self.dbName,'ten_elements_dups')
    d = conn.GetData(randomAccess=1,removeDups=-1)
    assert tuple(d[0])==(0,11)
    assert tuple(d[2])==(2,21)
    assert len(d)==20
    self.failUnlessRaises(IndexError,lambda:d[21])

    # repeat that test to make sure the table argument works
    conn = DbConnect(self.dbName,'ten_elements')
    d = conn.GetData(table='ten_elements_dups',randomAccess=1,removeDups=-1)
    assert tuple(d[0])==(0,11)
    assert tuple(d[2])==(2,21)
    assert len(d)==20
    self.failUnlessRaises(IndexError,lambda:d[21])

  def testGetData4(self):
    """ non random access
    """
    conn = DbConnect(self.dbName,'ten_elements')
    d = conn.GetData(randomAccess=0)
    self.failUnlessRaises(TypeError,lambda : len(d))
      
    rs = []
    for thing in d:
      rs.append(thing)
    assert len(rs)==10
    assert tuple(rs[0])==(0,11)
    assert tuple(rs[2])==(4,31)
    
  def testGetData5(self):
    """ using a RandomAccessDbResultSet with a Transform
    """
    fn = lambda x:(x[0],x[1]*2)
    conn = DbConnect(self.dbName,'ten_elements')
    d = conn.GetData(randomAccess=1,transform=fn)

    assert tuple(d[0])==(0,22),str(d[0])
    assert tuple(d[2])==(4,62)
    assert len(d)==10
    self.failUnlessRaises(IndexError,lambda:d[11])
    
  def testGetData6(self):
    """ using a DbResultSet with a Transform
    """
    fn = lambda x:(x[0],x[1]*2)
    conn = DbConnect(self.dbName,'ten_elements')
    d = conn.GetData(randomAccess=0,transform=fn)
    self.failUnlessRaises(TypeError,lambda : len(d))
    rs = []
    for thing in d:
      rs.append(thing)
    assert len(rs)==10
    assert tuple(rs[0])==(0,22)
    assert tuple(rs[2])==(4,62)
    

  def testInsertData(self):
    """ tests InsertData and InsertColumnData functionalities """
    newTblName = 'NEW_TABLE'
    conn = DbConnect(self.dbName)
    try:
      conn.GetCursor().execute('drop table %s'%(newTblName))
    except:
      pass
    conn.Commit()
    conn.AddTable(newTblName,'id int,val1 int, val2 int')
    for i in range(10):
      conn.InsertData(newTblName,(i,i+1,2*i))
    conn.Commit()
    d = conn.GetData(table=newTblName)
    assert len(d)==10

    d = None
    try:
      conn.GetCursor().execute('drop table %s'%(newTblName))
    except:
      assert 0,'drop table failed'
    
    
    
if __name__ == '__main__':
  unittest.main()

