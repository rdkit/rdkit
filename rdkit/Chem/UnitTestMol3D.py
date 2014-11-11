# $Id$
#
"""unit testing code for 3D stuff

"""
from rdkit import RDConfig
import unittest,os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import TorsionFingerprints

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
    self.assertAlmostEqual(rms2, 0.0, 4)


  def testConformerRMSMatrix(self):
    m1 = Chem.MolFromSmiles('CNc(n2)nc(C)cc2Nc(cc34)ccc3[nH]nc4')
    cids = AllChem.EmbedMultipleConfs(m1,3)

    m2 = Chem.MolFromSmiles('CNc(n2)nc(C)cc2Nc(cc34)ccc3[nH]nc4')
    m2.AddConformer(m1.GetConformer(id=0))

    # test that the RMS matrix has the correct size
    rmat = AllChem.GetConformerRMSMatrix(m1)
    self.assertEqual(len(rmat), 3)

    # test that the elements are in the right order
    self.assertAlmostEqual(rmat[0], AllChem.GetBestRMS(m1, m2, 1, 0), 3)
    self.assertAlmostEqual(rmat[1], AllChem.GetBestRMS(m1, m2, 2, 0), 3)

    # test the prealigned option
    rmat2 = AllChem.GetConformerRMSMatrix(m1, prealigned=True)
    self.assertAlmostEqual(rmat[0], rmat2[0])

  def testTorsionFingerprints(self):
    # we use the xray structure from the paper (JCIM, 52, 1499, 2012): 1DWD
    refFile = os.path.join(RDConfig.RDCodeDir,'Chem','test_data','1DWD_ligand.pdb')
    ref = Chem.MolFromSmiles('NC(=[NH2+])c1ccc(C[C@@H](NC(=O)CNS(=O)(=O)c2ccc3ccccc3c2)C(=O)N2CCCCC2)cc1')
    mol = Chem.MolFromPDBFile(refFile)
    mol = AllChem.AssignBondOrdersFromTemplate(ref, mol)
    
    # the torsion lists
    tors_list, tors_list_rings = TorsionFingerprints.CalculateTorsionLists(mol)
    self.assertEqual(len(tors_list), 11)
    self.assertEqual(len(tors_list_rings), 4)
    self.assertAlmostEqual(tors_list[-1][1], 180.0, 4)
    tors_list, tors_list_rings = TorsionFingerprints.CalculateTorsionLists(mol, maxDev='spec')
    self.assertAlmostEqual(tors_list[-1][1], 90.0, 4)
    self.assertRaises(ValueError, TorsionFingerprints.CalculateTorsionLists, mol, maxDev='test')
    tors_list, tors_list_rings = TorsionFingerprints.CalculateTorsionLists(mol, symmRadius=0)
    self.assertEqual(len(tors_list[0][0]), 2)

    # the weights
    weights = TorsionFingerprints.CalculateTorsionWeights(mol)
    self.assertAlmostEqual(weights[4], 1.0)
    self.assertEqual(len(weights),len(tors_list+tors_list_rings))
    weights = TorsionFingerprints.CalculateTorsionWeights(mol, 15, 14)
    self.assertAlmostEqual(weights[3], 1.0)
    self.assertRaises(ValueError, TorsionFingerprints.CalculateTorsionWeights, mol, 15, 3)

    # the torsion angles
    tors_list, tors_list_rings = TorsionFingerprints.CalculateTorsionLists(mol)
    torsions = TorsionFingerprints.CalculateTorsionAngles(mol, tors_list, tors_list_rings)
    self.assertEqual(len(weights), len(torsions))
    self.assertAlmostEqual(torsions[2][0], 232.5346, 4)

    # the torsion fingerprint deviation
    tfd = TorsionFingerprints.CalculateTFD(torsions, torsions)
    self.assertAlmostEqual(tfd, 0.0)
    refFile = os.path.join(RDConfig.RDCodeDir,'Chem','test_data','1PPC_ligand.pdb')
    mol2 = Chem.MolFromPDBFile(refFile)
    mol2 = AllChem.AssignBondOrdersFromTemplate(ref, mol2)
    torsions2 = TorsionFingerprints.CalculateTorsionAngles(mol2, tors_list, tors_list_rings)
    weights = TorsionFingerprints.CalculateTorsionWeights(mol)
    tfd = TorsionFingerprints.CalculateTFD(torsions, torsions2, weights=weights)
    self.assertAlmostEqual(tfd, 0.0645,4)
    tfd = TorsionFingerprints.CalculateTFD(torsions, torsions2)
    self.assertAlmostEqual(tfd, 0.1680,4)

if __name__ == '__main__':
  unittest.main()

