#
#  Copyright (C) 2003-2006  Greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""unit testing code for fingerprinting

"""
import unittest

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols

class TestCase(unittest.TestCase):

  def test1(self):
    # FIX: test HashAtom
    pass

  def test2(self):
    # FIX: test HashBond
    pass

  def test3(self):
    # FIX: test HashPath
    pass

  def test4(self):
    """ check containing mols, no Hs, no valence """
    tgts = [('CCC(O)C(=O)O', ('CCC', 'OCC', 'OCC=O', 'OCCO', 'CCCC', 'OC=O', 'CC(O)C')), ]
    for smi, matches in tgts:
      m = Chem.MolFromSmiles(smi)
      fp1 = Chem.RDKFingerprint(m, 2, 7, 9192, 4, 0)
      _ = fp1.GetOnBits()
      for match in matches:
        m2 = Chem.MolFromSmiles(match)
        fp2 = Chem.RDKFingerprint(m2, 2, 7, 9192, 4, 0)
        v1, _ = DataStructs.OnBitProjSimilarity(fp2, fp1)
        self.assertAlmostEqual(v1, 1, 'substruct %s not properly contained in %s' % (match, smi))

  def test5(self):
    """ check containing mols, use Hs, no valence """
    tgts = [('CCC(O)C(=O)O', ('O[CH-][CH2-]', 'O[CH-][C-]=O')), ]
    for smi, matches in tgts:
      m = Chem.MolFromSmiles(smi)
      fp1 = Chem.RDKFingerprint(m, 2, 7, 9192, 4, 1)
      _ = fp1.GetOnBits()
      for match in matches:
        m2 = Chem.MolFromSmiles(match)
        fp2 = Chem.RDKFingerprint(m2, 2, 7, 9192, 4, 1)
        v1, _ = DataStructs.OnBitProjSimilarity(fp2, fp1)
        self.assertAlmostEqual(v1, 1, 'substruct %s not properly contained in %s' % (match, smi))

  def test6(self):
    """ check that the bits in a signature of size N which has been folded in half
      are the same as those in a signature of size N/2 """
    smis = ['CCC(O)C(=O)O', 'c1ccccc1', 'C1CCCCC1', 'C1NCCCC1', 'CNCNCNC']
    for smi in smis:
      m = Chem.MolFromSmiles(smi)
      fp1 = Chem.RDKFingerprint(m, 2, 7, 4096)
      fp2 = DataStructs.FoldFingerprint(fp1, 2)
      fp3 = Chem.RDKFingerprint(m, 2, 7, 2048)
      self.assertEqual(tuple(fp2.GetOnBits()), tuple(fp3.GetOnBits()))

      fp2 = DataStructs.FoldFingerprint(fp2, 2)
      fp3 = Chem.RDKFingerprint(m, 2, 7, 1024)
      self.assertEqual(tuple(fp2.GetOnBits()), tuple(fp3.GetOnBits()))

      fp2 = DataStructs.FoldFingerprint(fp1, 4)
      self.assertEqual(tuple(fp2.GetOnBits()), tuple(fp3.GetOnBits()))

  def testGithub1747(self):
    """ test github #1747: deprecated apply() function causes GetRDKFingerprint
        to fail in Python 3 """
    fp = FingerprintMols.GetRDKFingerprint(Chem.MolFromSmiles('CCO'))
    self.assertNotEqual(0,fp.GetNumOnBits())

if __name__ == '__main__':  # pragma: nocover
  unittest.main()
