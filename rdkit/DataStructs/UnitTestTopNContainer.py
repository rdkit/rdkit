# $Id$
#
# Copyright (C) 2003-2006 Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
import random
import unittest

from io import StringIO

from rdkit.DataStructs.TopNContainer import TopNContainer, _exampleCode
from rdkit.TestRunner import redirect_stdout


class TestCase(unittest.TestCase):

    def test1(self):
        # simple test with a known answer
        cont = TopNContainer(4)
        for foo in range(10):
            cont.Insert(foo, str(foo))
        assert cont.GetPts() == list(range(6, 10))
        assert cont.GetExtras() == [str(x) for x in range(6, 10)]

    def test2(self):
        # larger scale random test
        cont = TopNContainer(50)
        for _ in range(1000):
            cont.Insert(random.random())
        vs = cont.GetPts()
        last = vs.pop(0)
        while vs:
            assert vs[0] >= last
            last = vs.pop(0)

    def test3(self):
        # random test with extras
        cont = TopNContainer(10)
        for _ in range(100):
            v = random.random()
            cont.Insert(v, v + 1)
        vs = cont.GetExtras()
        last = vs.pop(0)
        while vs:
            assert vs[0] >= last
            last = vs.pop(0)

    def test4(self):
        # random test with extras and getitem
        cont = TopNContainer(10)
        for i in range(100):
            v = random.random()
            cont.Insert(v, v + 1)
        lastV, lastE = cont[0]
        for i in range(1, len(cont)):
            v, e = cont[i]
            assert v >= lastV
            assert e >= lastE
            lastV, lastE = v, e

    def test5(self):
        # random test with extras and getitem, include reverse
        cont = TopNContainer(10)
        for i in range(100):
            v = random.random()
            cont.Insert(v, v + 1)
        cont.reverse()
        lastV, lastE = cont[0]
        for i in range(1, len(cont)):
            v, e = cont[i]
            assert v <= lastV
            assert e <= lastE
            lastV, lastE = v, e

    def test_keepAll(self):
        # simple test with a known answer where we keep all
        cont = TopNContainer(-1)
        for i in range(10):
            cont.Insert(9 - i, str(9 - i))
            self.assertEqual(len(cont), i + 1)
        assert cont.GetPts() == list(range(10))
        assert cont.GetExtras() == [str(x) for x in range(10)]

    def test_exampleCode(self):
        # We make sure that the example code runs
        f = StringIO()
        with redirect_stdout(f):
            _exampleCode()
        s = f.getvalue()
        self.assertIn('[58, 75, 78, 84]', s)


if __name__ == '__main__':  # pragma: nocover
    unittest.main()
