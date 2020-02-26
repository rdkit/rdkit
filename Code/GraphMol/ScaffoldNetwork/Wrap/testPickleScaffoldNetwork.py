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

from rdkit import Chem
from rdkit import RDConfig
from rdkit import rdBase
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


if __name__ == '__main__':
  unittest.main()
