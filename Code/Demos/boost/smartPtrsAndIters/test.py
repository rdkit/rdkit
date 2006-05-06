# Copyright Rational Discovery LLC 2005
# Distributed under the Boost Software License, Version 1.0. (See
# accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)
import SPtrTestModule as TestModule
import unittest


class TestCase(unittest.TestCase):
  def setUp(self):
    pass
  def test1(self):
    obj = TestModule.DemoKlass(3);
    self.failUnless(obj.GetVal()==3)
  def test2(self):
    obj = TestModule.buildPtr(3);
    self.failUnless(obj.GetVal()==3)
  def test3(self):
    obj = TestModule.buildSPtr(3);
    self.failUnless(obj.GetVal()==3)
  def test4(self):
    sz = 5
    vect = TestModule.buildPtrVector(sz)
    self.failUnless(len(vect)==sz)
    for i in range(sz):
      self.failUnless(vect[i].GetVal()==i)
  def test5(self):
    sz = 5
    vect = TestModule.buildSPtrVector(sz)
    self.failUnless(len(vect)==sz)
    for i in range(sz):
      self.failUnless(vect[i].GetVal()==i)
  def test5b(self):
    sz = 5
    vect = TestModule.buildSPtrVector(sz)
    self.failUnless(len(vect)==sz)
    p = 0
    for itm in vect:
      self.failUnless(itm.GetVal()==p)
      p+=1
  def test6(self):
    sz = 5
    cont = TestModule.DemoContainer(sz)
    p = 0
    for itm in cont:
      self.failUnless(itm.GetVal()==p)
      p+=1
    self.failUnless(p==sz)
if __name__ == '__main__':
  unittest.main()

