#
# Copyright (C) 2019 Greg Landrum and T5 Informatics GmbH
#  All Rights Reserved
#
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#

import os
import pickle
import unittest

from rdkit import Chem, RDConfig, rdBase
from rdkit.Chem.Scaffolds import rdScaffoldNetwork

rdBase.DisableLog("rdApp.info")


class TestNbScaffoldNetworkPickle(unittest.TestCase):
  """
  Tests for the nanobind ScaffoldNetwork serialization format.

  RDKit must have been built with Boost serialization support
  for these tests to pass.
  """

  def _makeNetwork(self):
    smis = ["c1ccccc1CC1NC(=O)CCC1", "c1cccnc1CC1NC(=O)CCC1"]
    ms = [Chem.MolFromSmiles(x) for x in smis]
    params = rdScaffoldNetwork.ScaffoldNetworkParams()
    params.includeScaffoldsWithoutAttachments = False
    return rdScaffoldNetwork.CreateScaffoldNetwork(ms, params)

  def test2PickleStateIsBytes(self):
    """__getstate__ must return a tuple containing bytes, not str.

    The nanobind wrap uses nb::bytes so that binary serialization formats
    remain safe. This test confirms the contract.
    """
    net = self._makeNetwork()
    state = net.__getstate__()
    self.assertIsInstance(state, tuple)
    self.assertEqual(len(state), 1)
    self.assertIsInstance(state[0], bytes)

  def test3LegacyBoostPickle(self):
    """The nanobind wrap must load pickles created by the boost.python wrap.

    scaffold_network_boost.pkl was generated with the boost.python
    rdScaffoldNetwork from the same two SMILES used in _makeNetwork().
    """
    pkl_path = os.path.join(RDConfig.RDBaseDir,
                            'Code/GraphMol/ScaffoldNetwork/Wrap/testData',
                            'scaffold_network_boost.pkl')
    with open(pkl_path, 'rb') as f:
      net2 = pickle.load(f)

    net = self._makeNetwork()
    self.assertEqual(len(net2.nodes), len(net.nodes))
    self.assertEqual(sorted(net2.nodes), sorted(net.nodes))
    self.assertEqual(len(net2.edges), len(net.edges))


if __name__ == '__main__':
  unittest.main()
