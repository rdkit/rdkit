# $Id$
#
#  Copyright (C) 2001-2006  greg Landrum and Rational Discovery LLC
#
#   @@ All Rights Reserved @@
#  This file is part of the RDKit.
#  The contents are covered by the terms of the BSD license
#  which is included in the file license.txt, found at the root
#  of the RDKit source tree.
#
"""testing code inspired by old bugs

The bugs were in the OELib code, so these are maybe no longer
relevant... but tests are tests

"""
from rdkit import RDConfig
import unittest, os
import pickle
from rdkit import Chem
from rdkit.Chem import AllChem


def feq(n1, n2, tol=1e-4):
  return abs(n1 - n2) <= tol


class TestCase(unittest.TestCase):

  def testBug12a(self):
    from rdkit.Chem import MolSurf
    inD = [
      ('OC(=O)[CH](CC1=CC=CC=C1)C2=CC=CC=C2', 37.3),
      ('OC(=O)C(C1=CC=CC=C1)(C2=CC=CC=C2)C3=CC=CC=C3', 37.3),
      ('CCC(CC)(CC)[CH](OC(=O)C1=C(C=CC=C1)C(O)=O)C2=CC=CC=C2', 63.6),
      ('C[C](O)([CH](C(O)=O)C1=CC=CC=C1)C2=CC=CC=C2', 57.53),
      ('C[CH]([CH](C(O)=O)C1=CC=CC=C1)C2=CC=CC=C2', 37.3),
      ('OC(=O)CBr', 37.3),
      ('OC(=O)CCl', 37.3),
      ('OC(=O)C=CC(=O)C1=CC=CC=C1', 54.37),
      ('NC1=C(C=CC=C1)C(O)=O', 63.32),
      ('OC(=O)C1=CC=CC=C1', 37.3),
      ('CN(C)C(=N)NC1=NC(=C2C=C(Cl)C=CC2=N1)C.O[N+]([O-])=O', 128.27),
      ('CCN(CC)C(=N)NC1=NC(=C2C=C(Cl)C=CC2=N1)C.O[N+]([O-])=O', 128.27),
      ('ON(O)NC(=N)NN=C1C(=O)NC2=C1C=CC=C2', 133.07),
      ('NC1=CC=C(C=C1)C=NNC(=N)NN(O)O', 129.99),
      ('CC(=O)NC1=CC=C(C=C1)C=NNC(=N)NN(O)O', 133.07),
      ('COC1=CC=C(C=C1)C=NNC(=N)NN(O)O', 113.2),
      ('ON(O)NC(=N)NN=CC1=CC=CC=C1', 103.97),
      ('ON(O)NC(=N)NN=CC=CC1=CC=CC=C1', 103.97),
      ('ON(O)NC(=N)NN=CC1=C(Cl)C=C(Cl)C=C1', 103.97),
      ('CC(C)=CCCC(C)=CC=NNC(=N)NN(O)O', 103.97),
      ('CN(C)C1=CC=C(C=C1)C=NNC(=N)NN(O)O', 107.21),
      ('ON(O)NC(=N)NN=CC1=CC=CO1', 117.11),
      ('ON(O)NC(=N)NN=CC1=CC=C(O)C=C1', 124.2),
      ('CC(C)C1=CC=C(C=C1)C=NNC(=N)NN(O)O', 103.97),
      ('COC1=C(C=CC=C1)C=NNC(=N)NN(O)O', 113.2),
      ('ON(O)NC(=N)NN=CC1=C(C=CC=C1)[N+]([O-])=O', 147.11),
      ('ON(O)NC(=N)NN=CC1=CC=C(C=C1)[N+]([O-])=O', 147.11),
      ('ON(O)NC(=N)NN=CC1=C(O)C=CC(=C1)[N+]([O-])=O', 167.34),
      ('ON(O)NC(=N)NN=CC1=CC=NC=C1', 116.86),
      ('ON(O)NC(=N)NN=CC1=CC=CC=N1', 116.86),
      ('ON(O)NC(=N)NN=CC1=CC=CN=C1', 116.86),
    ]
    for smi, val in inD:
      mol = Chem.MolFromSmiles(smi)
      v = MolSurf.TPSA(mol)
      assert feq(v, val), 'bad TPSA (%f != %f) for smiles: %s' % (v, val, smi)

  def testBug12b(self):
    """ failures for Bug12 which are actually related to Bug14

    """
    from rdkit.Chem import MolSurf
    inD = [('[O-][N+](=O)C1=CNC(=N)S1', 82.78), ]
    for smi, val in inD:
      mol = Chem.MolFromSmiles(smi)
      v = MolSurf.TPSA(mol)
      assert feq(v, val), 'bad TPSA (%f != %f) for smiles: %s' % (v, val, smi)

  def testBug14(self):
    """ 
    
    """
    smi = '[O-][N+](=O)C1=CNC(=N)S1'
    mol = Chem.MolFromSmiles(smi)
    at = mol.GetAtomWithIdx(5)
    assert at.GetHybridization() == Chem.HybridizationType.SP2, 'bad hyb'
    assert at.GetTotalNumHs() == 1, 'bad H count'

    mol = Chem.MolFromSmiles(smi)
    at = mol.GetAtomWithIdx(5)
    assert at.GetTotalNumHs() == 1, 'bad H count'
    assert at.GetHybridization() == Chem.HybridizationType.SP2, 'bad hyb'

  def testGithub112(self):
    """
    problems with AllChem.GetBestRMS() and molecules with Hs
    
    """
    m0 = Chem.MolFromMolFile(
      os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data', 'github112_tgt.mol'), removeHs=False)
    m1 = Chem.MolFromMolFile(
      os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data', 'github112_qry.mol'), removeHs=False)
    rms = AllChem.GetBestRMS(m0, m1)
    self.assertAlmostEqual(rms, 0.456, 3)


if __name__ == '__main__':
  unittest.main()
