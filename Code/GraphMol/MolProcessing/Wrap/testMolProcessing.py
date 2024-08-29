#
#  Copyright (C) 2024 Greg Landrum and other RDKit contributors
#   @@ All Rights Reserved @@
#
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.

import unittest

#
from rdkit import Chem
from rdkit.Chem import rdMolProcessing
from rdkit import RDConfig


class TestCase(unittest.TestCase):

  def setUp(self):
    self.smiFile = RDConfig.RDBaseDir + '/Regress/Data/zinc.leads.500.q.smi'
    self.sdFile = RDConfig.RDBaseDir + "/Data/NCI/first_200.props.sdf"

  def test1(self):
    fps = rdMolProcessing.GetFingerprintsForMolsInFile(self.smiFile)
    self.assertEqual(len(fps), 499)
    fps = rdMolProcessing.GetFingerprintsForMolsInFile(self.sdFile)
    self.assertEqual(len(fps), 200)


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
