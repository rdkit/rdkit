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
from rdkit.Chem import rdFingerprintGenerator
from rdkit import DataStructs
from rdkit import RDConfig


class TestCase(unittest.TestCase):

  def setUp(self):
    self.smiFile = RDConfig.RDBaseDir + '/Regress/Data/zinc.leads.500.q.smi'
    self.sdFile = RDConfig.RDBaseDir + "/Data/NCI/first_200.props.sdf"

  def test1(self):
    fpg = rdFingerprintGenerator.GetMorganGenerator()
    fps = rdMolProcessing.GetFingerprintsForMolsInFile(self.smiFile)
    self.assertEqual(len(fps), 499)
    with Chem.SmilesMolSupplier(self.smiFile, delimiter='\t') as suppl:
      mols = [next(suppl) for _ in range(3)]
    nfps = [fpg.GetFingerprint(m) for m in mols]
    self.assertEqual(DataStructs.TanimotoSimilarity(fps[0], fps[1]),
                     DataStructs.TanimotoSimilarity(nfps[0], nfps[1]))

    fps = rdMolProcessing.GetFingerprintsForMolsInFile(self.sdFile)
    self.assertEqual(len(fps), 200)
    with Chem.SDMolSupplier(self.sdFile) as suppl:
      mols = [next(suppl) for _ in range(3)]
    nfps = [fpg.GetFingerprint(m) for m in mols]
    self.assertAlmostEqual(DataStructs.TanimotoSimilarity(fps[0], fps[1]), 0.0638, places=3)

  def test2(self):
    fpg = rdFingerprintGenerator.GetMorganGenerator(radius=2)

    fps = rdMolProcessing.GetFingerprintsForMolsInFile(self.smiFile, generator=fpg)
    self.assertEqual(len(fps), 499)
    with Chem.SmilesMolSupplier(self.smiFile, delimiter='\t') as suppl:
      mols = [next(suppl) for _ in range(3)]
    nfps = [fpg.GetFingerprint(m) for m in mols]
    self.assertEqual(DataStructs.TanimotoSimilarity(fps[0], fps[1]),
                     DataStructs.TanimotoSimilarity(nfps[0], nfps[1]))


if __name__ == '__main__':  # pragma: nocover
  unittest.main()
