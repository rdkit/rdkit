# $Id$
#
#  Copyright (C) 2002-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved  @@
#
"""unit testing code for the database utilities

"""
import RDConfig
import unittest,os
from Dbase import DbUtils
from Dbase.DbConnection import DbConnect

class TestCase(unittest.TestCase):
  def setUp(self):
    #print '\n%s: '%self.shortDescription(),
    self.baseDir = os.path.join(RDConfig.RDCodeDir,'Dbase','testData')
    if not RDConfig.usePgSQL:
      self.dbName = os.path.join(self.baseDir,'TEST.GDB')
    else:
      self.dbName = "::RDTests"
    self.colHeads=('int_col','floatCol','strCol')
    self.colTypes=('integer','float','string')

  def _confirm(self,tblName):
    conn = DbConnect(self.dbName,tblName)
    res = conn.GetColumnNamesAndTypes()
    assert len(res)==len(self.colHeads),'bad number of columns'
    names = [x[0] for x in res]
    for i in range(len(names)):
      assert names[i].upper()==self.colHeads[i].upper(),'bad column head'
    types = [x[1] for x in res]
    for i in range(len(types)):
      assert types[i]==self.colTypes[i],'bad column type'
    
  def test1Txt(self):
    """ test reading from a text file """
    inF = open(os.path.join(self.baseDir,'dbtest.csv'),'r')
    tblName = 'fromtext'
    DbUtils.TextFileToDatabase(self.dbName,tblName,inF)
    self._confirm(tblName)

  def test2Excel(self):
    """ test the Excel functionality """
    import sys
    if sys.platform!='win32':
      return
    from Excel import ExcelWrapper
    w = ExcelWrapper.ExcelWrapper()
    xlFile = os.path.join(self.baseDir,'dbtest.xls')
    w.Workbooks.Open(xlFile)
    tblName = 'fromxl'
    DbUtils.ExcelToDatabase(self.dbName,tblName,wrapper=w,force=1)
    w.ActiveWorkbook.Close()
    if len(w.Workbooks)==0: w.Quit()
    self._confirm(tblName)

  def test3Txt(self):
    """ test reading from a text file including null markers"""
    inF = open(os.path.join(self.baseDir,'dbtest.nulls.csv'),'r')
    tblName = 'fromtext2'
    DbUtils.TextFileToDatabase(self.dbName,tblName,inF,nullMarker='NA')
    self._confirm(tblName)

  def test4Excel(self):
    """ test the Excel functionality including nulls """
    import sys
    if sys.platform!='win32':
      return
    from Excel import ExcelWrapper
    w = ExcelWrapper.ExcelWrapper()
    xlFile = os.path.join(self.baseDir,'dbtest.nulls.xls')
    w.Workbooks.Open(xlFile)
    tblName = 'fromxl2'
    DbUtils.ExcelToDatabase(self.dbName,tblName,wrapper=w,nullMarker='NA',force=1)
    w.ActiveWorkbook.Close()
    if len(w.Workbooks)==0: w.Quit()
    self._confirm(tblName)


  def testGetData1(self):
    """ basic functionality
    """
    d = DbUtils.GetData(self.dbName,'ten_elements',forceList=1)
    assert len(d)==10
    assert tuple(d[0])==(0,11)
    assert tuple(d[2])==(4,31)
    try:
      d[11]
    except IndexError:
      pass
    except:
      assert 0,'bad exception type raised'
    else:
      assert 0,'failed to raise expected exception'

  def testGetData2(self):
    """ using a RandomAccessDbResultSet
    """
    d = DbUtils.GetData(self.dbName,'ten_elements',forceList=0,randomAccess=1)
    assert tuple(d[0])==(0,11)
    assert tuple(d[2])==(4,31)
    assert len(d)==10
    try:
      d[11]
    except IndexError:
      pass
    except:
      assert 0,'bad exception type raised'
    else:
      assert 0,'failed to raise expected exception'

  def testGetData3(self):
    """ using a DbResultSet
    """
    d = DbUtils.GetData(self.dbName,'ten_elements',forceList=0,randomAccess=0)
    try:
      len(d)
    except TypeError:
      pass
    except:
      assert 0,'bad exception type raised'
    else:
      assert 0,'failed to raise expected exception'
    rs = []
    for thing in d:
      rs.append(thing)
    assert len(rs)==10
    assert tuple(rs[0])==(0,11)
    assert tuple(rs[2])==(4,31)
    
  def testGetData4(self):
    """ using a RandomAccessDbResultSet with a Transform
    """
    fn = lambda x:(x[0],x[1]*2)
    d = DbUtils.GetData(self.dbName,'ten_elements',forceList=0,randomAccess=1,
                        transform=fn)
    assert tuple(d[0])==(0,22)
    assert tuple(d[2])==(4,62)
    assert len(d)==10
    try:
      d[11]
    except IndexError:
      pass
    except:
      assert 0,'bad exception type raised'
    else:
      assert 0,'failed to raise expected exception'

    
  def testGetData5(self):
    """ using a DbResultSet with a Transform
    """
    fn = lambda x:(x[0],x[1]*2)
    d = DbUtils.GetData(self.dbName,'ten_elements',forceList=0,randomAccess=0,
                        transform=fn)
    try:
      len(d)
    except TypeError:
      pass
    except:
      assert 0,'bad exception type raised'
    else:
      assert 0,'failed to raise expected exception'
    rs = []
    for thing in d:
      rs.append(thing)
    assert len(rs)==10
    assert tuple(rs[0])==(0,22)
    assert tuple(rs[2])==(4,62)
    

    
    
    
if __name__ == '__main__':
  unittest.main()

