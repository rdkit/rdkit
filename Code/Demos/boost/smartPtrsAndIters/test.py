# Copyright Rational Discovery LLC 2005
# Distributed under the Boost Software License, Version 1.0. (See
# accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)
import unittest

import SPtrTestModule as TestModule


class TestCase(unittest.TestCase):

  def setUp(self):
    pass

  def test1(self):
    obj = TestModule.DemoKlass(3)
    self.assertTrue(obj.GetVal() == 3)

  def test2(self):
    obj = TestModule.buildPtr(3)
    self.assertTrue(obj.GetVal() == 3)

  def test3(self):
    obj = TestModule.buildSPtr(3)
    self.assertTrue(obj.GetVal() == 3)

  def test4(self):
    sz = 5
    vect = TestModule.buildPtrVector(sz)
    self.assertTrue(len(vect) == sz)
    for i in range(sz):
      self.assertTrue(vect[i].GetVal() == i)

  def test5(self):
    sz = 5
    vect = TestModule.buildSPtrVector(sz)
    self.assertTrue(len(vect) == sz)
    for i in range(sz):
      self.assertTrue(vect[i].GetVal() == i)

  def test5b(self):
    sz = 5
    vect = TestModule.buildSPtrVector(sz)
    self.assertTrue(len(vect) == sz)
    p = 0
    for itm in vect:
      self.assertTrue(itm.GetVal() == p)
      p += 1

  def test6(self):
    sz = 5
    cont = TestModule.DemoContainer(sz)
    p = 0
    for itm in cont:
      self.assertTrue(itm.GetVal() == p)
      p += 1
    self.assertTrue(p == sz)


if __name__ == '__main__':
  unittest.main()
