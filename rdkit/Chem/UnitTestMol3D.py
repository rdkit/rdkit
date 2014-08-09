# $Id$
#
"""unit testing code for 3D stuff

"""
from rdkit import RDConfig
import unittest,os
from rdkit import Chem
from rdkit.Chem import AllChem

def feq(n1,n2,tol=1e-4):
  return abs(n1-n2)<=tol

class TestCase(unittest.TestCase):

  def testConformerRMS(self):
    m1 = Chem.MolFromSmiles('CNc(n2)nc(C)cc2Nc(cc34)ccc3[nH]nc4')
    cids = AllChem.EmbedMultipleConfs(m1,2)

    m2 = Chem.MolFromSmiles('CNc(n2)nc(C)cc2Nc(cc34)ccc3[nH]nc4')
    m2.AddConformer(m1.GetConformer(id=1))

    # test that the prealigned flag is working
    rms1 = AllChem.GetConformerRMS(m1, 0, 1, prealigned=True)
    rms2 = AllChem.GetConformerRMS(m1, 0, 1, prealigned=False)
    self.assertTrue((rms1>rms2))

    # test that RMS is the same as calculated by AlignMol()
    self.assertAlmostEqual(rms2, AllChem.GetBestRMS(m2, m1, 1, 0), 3)

    # the RMS with itself must be zero
    rms2 = AllChem.GetConformerRMS(m1, 0, 0, prealigned=True)
    self.assertTrue(feq(rms2, 0.0))


  def testConformerRMSMatrix(self):
    m1 = Chem.MolFromSmiles('CNc(n2)nc(C)cc2Nc(cc34)ccc3[nH]nc4')
    cids = AllChem.EmbedMultipleConfs(m1,3)

    m2 = Chem.MolFromSmiles('CNc(n2)nc(C)cc2Nc(cc34)ccc3[nH]nc4')
    m2.AddConformer(m1.GetConformer(id=0))

    # test that the RMS matrix has the correct size
    rmat = AllChem.GetConformerRMSMatrix(m1)
    self.assertTrue(len(rmat), 3)

    # test that the elements are in the right order
    self.assertAlmostEqual(rmat[0], AllChem.GetBestRMS(m1, m2, 1, 0), 3)
    self.assertAlmostEqual(rmat[1], AllChem.GetBestRMS(m1, m2, 2, 0), 3)


if __name__ == '__main__':
  unittest.main()

