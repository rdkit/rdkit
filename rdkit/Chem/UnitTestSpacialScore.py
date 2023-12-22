import unittest

from rdkit import Chem
from rdkit.Chem.SpacialScore import SPS

# The tests reproduce the results from the test in the original repository (entry 0) 
# (https://github.com/frog2000/Spacial-Score/blob/main/test_spacial_score.py)
# as well as the entries from Table 1 of the manuscript (entries 1-)
# https://doi.org/10.1021/acs.jmedchem.3c00689
#
# SMILES, expected SPS, expected nSPS:
TEST_CASES = [
  (r"C/C=C\C1C=CCC(Br)C1C2=CC(C#C)=CC=C2", 491, 27.28),
  ("C=CCCBr", 37, 7.40),
  ("CCCCC", 42, 8.40),
  ("CCCCBr", 42, 8.40),
  ("CC(C)CBr", 48, 9.60),
  ("C/C=C/CBr", 50, 10.00),
  ("BrC(C)(C)(C)", 60, 12.00),
  ("CC[C@@H](C)Br", 75, 15.00),
  ("c1ccccc1", 48, 8.00),
  ("c1cnccc1", 48, 8.00),
  ("C1=CC=CCC1", 112, 18.67),
  ("C1=CCCCC1", 128, 21.33),
  ("C1CCCCC1", 144, 24.00),
  ("C1CCCC1C", 153, 25.50),
  ("C1CCC1(C)C", 174, 29.00),
  ("O=C1COCCN1c4ccc(N3C[C@H](CNC(=O)c2ccc(Cl)s2)OC3=O)cc4", 563, 19.41),
  ("CC(C)CCC[C@@H](C)[C@H]3CC[C@H]4C2CC=C1C[C@@H](O)CC[C@]1(C)[C@H]2CC[C@]34C", 1303, 46.54),
  
  ('CC=CC', 38, 9.50),
  ('C/C=C/C', 38, 9.50),
  ('C=C1CCCCC1', 158, 22.57),
  ('C=C1CCCNC1', 158, 22.57),
  ('CC=C1CCCCC1', 167, 20.88),
  ('C/C=C\\1CCCCC1', 167, 20.88),
  ('CC=C1CCCNC1', 211, 26.38),
  ('CC=C(C1CCCNC1)C2CCCNC2', 485, 32.33),
  ('C/C=C(C1CCCNC1)/C2CCCNC2', 485, 32.33),
  ('CC=C(C1CCCNC1)C2CCCCN2', 511, 34.07),
  ('C/C=C(C1CCCNC1)/C2CCCCN2', 511, 34.07),
  ('CC1=CCCCC1', 151, 21.57),
]


class TestCase(unittest.TestCase):

  def testVersion(self):
    self.assertEqual(SPS.version, '1.0.0',
                     msg='SpacialScore version has changed. Update the tests if required.')

  def testCases(self):
    for idx, tc in enumerate(TEST_CASES):
      mol = Chem.MolFromSmiles(tc[0])
      sps = SPS(mol, False)
      nsps = SPS(mol, True)
      self.assertEqual(sps, tc[1], msg=f"SPS {sps} not equal to expected value of {tc[1]} for entry {idx}")
      self.assertAlmostEqual(nsps, tc[2], places=2, msg=f"nSPS {nsps} not close to expected value of {tc[2]} for entry {idx}.")



if __name__ == '__main__':  # pragma: nocover
  unittest.main()
