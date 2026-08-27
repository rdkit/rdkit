#
# Copyright (C) 2019 Greg Landrum and T5 Informatics GmbH
#  All Rights Reserved
#
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#

import pickle
import unittest

from rdkit import Chem, rdBase
from rdkit.Chem.Scaffolds import rdScaffoldNetwork

rdBase.DisableLog("rdApp.info")


class TestScaffoldNetwork(unittest.TestCase):
  """
  RDKit must have been built with Boost serialization support
  for this test to pass.
  """

  def test1Pickle(self):
    smis = ["c1ccccc1CC1NC(=O)CCC1", "c1cccnc1CC1NC(=O)CCC1"]
    ms = [Chem.MolFromSmiles(x) for x in smis]
    params = rdScaffoldNetwork.ScaffoldNetworkParams()
    params.includeScaffoldsWithoutAttachments = False
    net = rdScaffoldNetwork.CreateScaffoldNetwork(ms, params)
    self.assertEqual(len(net.nodes), 7)
    self.assertEqual(len(net.edges), 7)

    pkl = pickle.dumps(net)
    net2 = pickle.loads(pkl)
    self.assertEqual(len(net2.nodes), 7)
    self.assertEqual(len(net2.edges), 7)
    self.assertEqual(list(net2.nodes), list(net.nodes))
    self.assertEqual([str(x) for x in net2.edges], [str(x) for x in net.edges])

    # we must be able to load pickles created by the boost.python wrappers.'''
    pkl = b'\x80\x04\x95N\x01\x00\x00\x00\x00\x00\x00\x8c&rdkit.Chem.Scaffolds.rdScaffoldNetwork\x94\x8c\x0fScaffoldNetwork\x94\x93\x94B\x01\x01\x00\x0022 serialization::archive 19 0 1 0 0 7 0 21 O=C1CCCC(Cc2ccccc2)N1 9 *c1ccccc1 15 **1:*:*:*:*:*:1 13 *C1CCCC(=O)N1 13 **1****(=*)*1 21 O=C1CCCC(Cc2ccccn2)N1 9 *c1ccccn1 7 0 1 1 2 2 2 1 1 7 0 1 1 2 2 2 1 1 0 0 7 0 0 0 0 1 1 1 2 2 0 3 1 3 4 2 5 6 1 6 2 2 5 3 1\x94\x85\x94R\x94}\x94\x85\x94b.'
    net2 = pickle.loads(pkl)
    self.assertEqual(len(net2.nodes), 7)
    self.assertEqual(len(net2.edges), 7)
    self.assertEqual(list(net2.nodes), list(net.nodes))
    self.assertEqual([str(x) for x in net2.edges], [str(x) for x in net.edges])


if __name__ == '__main__':
  unittest.main()
