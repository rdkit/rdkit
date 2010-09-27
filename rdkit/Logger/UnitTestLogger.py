# $Id$
#
#       Copyright (c) 2001-2006, Greg Landrum and Rational Discovery LLC,
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
""" unit testing code for the logger

"""
import unittest
import Logger
import re
class Foo:
  """ a simple class
  """
  def __init__(self,aVal):
    self.a = aVal
    
  def method1(self,a,b,c='foo'):
    return 'method1'

  def method2(self):
    return 'method2'

  def test1(self):
    return 'test1'

  def test2(self):
    return 'test2'


class TestCase(unittest.TestCase):
  def testBasic(self):
    l = Logger.Logger(Foo,7)
    try:
      l.method1(1,2)
      l.method1(1,2,c='grm')
      l.method2()
      l.method1(7,6,'pizza')
      l.b = 3
      l.method2()
    except:
      ok = 0
    else:
      ok = 1
    assert ok,'basic calls failed'

    assert l.b == 3, '__setattr__ failed'
    res = l._LoggerGetLog()
    assert len(res) == 6, 'length of log (%d) wrong'%(len(res))

  def testPlayback(self):
    l = Logger.Logger(Foo,7)
    try:
      l.method1(1,2)
      l.method1(1,2,c='grm')
      l.method2()
      l.method1(7,6,'pizza')
      l.b = 3
      l.method2()
    except:
      ok = 0
    else:
      ok = 1
    assert ok,'basic calls failed'

    f = Foo(7)
    replay = Logger.replay(l._LoggerGetLog(),f)
    assert replay==['method1','method1','method2','method1',3,'method2'],\
           'replay results (%s) wrong'%(str(replay))
    assert f.b == 3, '__setattr__ failed'
      
  def testFlush(self):
    l = Logger.Logger(Foo,7,loggerFlushCommand='method2')
    try:
      l.method1(1,2)
      l.method1(1,2,c='grm')
      l.method1(7,6,'pizza')
      l.b = 3
    except:
      ok = 0
    else:
      ok = 1
    assert ok,'basic calls failed'
    
    res = l._LoggerGetLog()
    assert len(res)==4, 'length of log (%d) wrong'%(len(res))

    l.method2()
    res = l._LoggerGetLog()
    assert len(res)==0, 'length of log (%d) wrong'%(len(res))
    
  def testIgnore(self):
    e1 = re.compile('test*')
    l = Logger.Logger(Foo,7,loggerIgnore=[e1,'method2'])
    try:
      l.method1(1,2)
      l.method2()
      l.test1()
      l.test2()
      l.method1(1,2,c='grm')
      l.method1(7,6,'pizza')
      l.b = 3
    except:
      ok = 0
    else:
      ok = 1
    assert ok,'basic calls failed'
    
    res = l._LoggerGetLog()
    assert len(res)==4, 'length of log (%d) wrong'%(len(res))


def TestSuite():
  suite = unittest.TestSuite()
  suite.addTest(TestCase('testBasic'))
  suite.addTest(TestCase('testPlayback'))
  suite.addTest(TestCase('testFlush'))
  suite.addTest(TestCase('testIgnore'))
  return suite


if __name__ == '__main__':
  suite = TestSuite()
  unittest.TextTestRunner().run(suite)

