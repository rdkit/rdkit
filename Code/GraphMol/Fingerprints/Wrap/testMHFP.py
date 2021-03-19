"""
2019, Daniel Probst, Reymond Group @ University of Bern
 @@ All Rights Reserved @@
This file is part of the RDKit.
The contents are covered by the terms of the BSD license
which is included in the file license.txt, found at the root
of the RDKit source tree.
"""

from rdkit import Chem
from rdkit.Chem import rdMHFPFingerprint
import unittest


class TestCase(unittest.TestCase):

  def setUp(self):
    pass

  def testMHFPFingerprint(self):
    s = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    t = "Cn1cnc2c1c(=O)[nH]c(=O)n2C"

    m = Chem.MolFromSmiles(s)
    enc = rdMHFPFingerprint.MHFPEncoder(128, 42)

    self.assertEqual(len(enc.CreateShinglingFromSmiles(s, rings=False)), 42)
    self.assertEqual(len(enc.CreateShinglingFromSmiles(s, min_radius=0)), 58)

    sh_a = enc.CreateShinglingFromSmiles(s)
    sh_b = enc.CreateShinglingFromMol(m)

    self.assertEqual(len(sh_a), 44)
    self.assertEqual(list(sh_a), list(sh_b))

    fp_a = enc.EncodeSmiles(s)
    fp_b = enc.EncodeMol(m)

    self.assertEqual(list(fp_a), list(fp_b))

    fp_c = enc.EncodeSmiles(t)
    dist = rdMHFPFingerprint.MHFPEncoder.Distance(fp_a, fp_c)
    self.assertEqual(dist, 0.4609375)


if __name__ == "__main__":
  unittest.main()
