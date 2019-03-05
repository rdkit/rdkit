# $Id$
#
#  Copyright (C) 2007  greg Landrum
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#


import doctest
import unittest

from rdkit.VLib import Node, Filter, Output, Supply, Transform
from io import StringIO


def load_tests(loader, tests, ignore):
    """ Add the Doctests from the module """
    tests.addTests(doctest.DocTestSuite(Filter, optionflags=doctest.ELLIPSIS))
    tests.addTests(doctest.DocTestSuite(Node, optionflags=doctest.ELLIPSIS))
    tests.addTests(doctest.DocTestSuite(Output, optionflags=doctest.ELLIPSIS))
    tests.addTests(doctest.DocTestSuite(Supply, optionflags=doctest.ELLIPSIS))
    tests.addTests(doctest.DocTestSuite(Transform, optionflags=doctest.ELLIPSIS))
    return tests


class Test_VLib(unittest.TestCase):

    def test_SupplyNode(self):
        supplier = Supply.SupplyNode()
        self.assertEqual(supplier._contents, [])

        supplier = Supply.SupplyNode(contents=[1, 2, 3])
        self.assertRaises(ValueError, supplier.AddParent, None)

    def test_FilterNode(self):
        filt = Filter.FilterNode(func=lambda a, b: a + b < 5)
        suppl1 = Supply.SupplyNode(contents=[1, 2, 3, 3])
        suppl2 = Supply.SupplyNode(contents=[1, 2, 3, 1])
        filt.AddParent(suppl1)
        filt.AddParent(suppl2)
        self.assertEqual([x for x in filt], [(1, 1), (2, 2), (3, 1)])
        filt.reset()
        self.assertEqual(filt.Negate(), False)
        filt.SetNegate(True)
        self.assertEqual(filt.Negate(), True)
        self.assertEqual([x for x in filt], [(3, 3), ])
        filt.Destroy()

    def test_OutputNode(self):
        supplier1 = Supply.SupplyNode(contents=[1, 2, 3])
        supplier2 = Supply.SupplyNode(contents=['a', 'b', 'c'])

        sio = StringIO()
        node = Output.OutputNode(dest=sio, strFunc=lambda x: '{0[0]}-{0[1]} '.format(x))
        node.AddParent(supplier1)
        node.AddParent(supplier2)
        result = list(s for s in node)
        self.assertEqual(result, [(1, 'a'), (2, 'b'), (3, 'c')])
        self.assertEqual(sio.getvalue(), '1-a 2-b 3-c ')

        sio = StringIO()
        node = Output.OutputNode(dest=sio)
        node.AddParent(supplier1)
        result = list(s for s in node)
        self.assertEqual(result, [1, 2, 3])
        self.assertEqual(sio.getvalue(), '123')

    def test_VLibNode(self):

        def setupNodes():
            p1 = Node.VLibNode()
            p2 = Node.VLibNode()
            c1 = Node.VLibNode()
            c2 = Node.VLibNode()
            p1.AddChild(c1)
            p2.AddChild(c1)
            p2.AddChild(c2)
            return p1, p2, c1, c2

        p1, p2, c1, c2 = setupNodes()
        # p1 -> c1
        # p2 -> c1
        # p2 -> c2
        self.assertEqual(len(c1.GetParents()), 2)
        self.assertEqual(len(c2.GetParents()), 1)
        self.assertEqual(len(p1.GetChildren()), 1)
        self.assertEqual(len(p2.GetChildren()), 2)

        p1.Destroy()
        # p1
        # p2 -> c1
        # p2 -> c2
        self.assertEqual(len(c1.GetParents()), 1)
        self.assertEqual(len(c2.GetParents()), 1)
        self.assertEqual(len(p1.GetChildren()), 0)
        self.assertEqual(len(p2.GetChildren()), 2)

        p1, p2, c1, c2 = setupNodes()
        p1.Destroy(propagateDown=True)
        # p1, c1
        # p2 -> c2
        self.assertEqual(len(c1.GetParents()), 0)
        self.assertEqual(len(c2.GetParents()), 1)
        self.assertEqual(len(p1.GetChildren()), 0)
        self.assertEqual(len(p2.GetChildren()), 1)

        p1, p2, c1, c2 = setupNodes()
        p1.Destroy(propagateUp=True)
        # p1
        # p2 -> c1
        # p2 -> c2
        self.assertEqual(len(c1.GetParents()), 1)
        self.assertEqual(len(c2.GetParents()), 1)
        self.assertEqual(len(p1.GetChildren()), 0)
        self.assertEqual(len(p2.GetChildren()), 2)

        p1, p2, c1, c2 = setupNodes()
        c1.Destroy(propagateUp=True)
        # p1, c1, p2, c2
        self.assertEqual(len(c1.GetParents()), 0)
        self.assertEqual(len(c2.GetParents()), 0)
        self.assertEqual(len(p1.GetChildren()), 0)
        self.assertEqual(len(p2.GetChildren()), 0)

        p1, p2, c1, c2 = setupNodes()
        p1.Destroy(propagateDown=True)
        # p1, c1
        # p2 -> c2
        self.assertEqual(len(c1.GetParents()), 0)
        self.assertEqual(len(c2.GetParents()), 1)
        self.assertEqual(len(p1.GetChildren()), 0)
        self.assertEqual(len(p2.GetChildren()), 1)
        p1.Destroy(propagateDown=True)
        # p1, c1
        # p2 -> c2
        self.assertEqual(len(c1.GetParents()), 0)
        self.assertEqual(len(c2.GetParents()), 1)
        self.assertEqual(len(p1.GetChildren()), 0)
        self.assertEqual(len(p2.GetChildren()), 1)


if __name__ == '__main__':  # pragma: nocover
    unittest.main()
