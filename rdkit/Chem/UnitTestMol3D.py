"""unit testing code for 3D stuff

"""
import os
import random
import unittest

from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, TorsionFingerprints

class TestCase(unittest.TestCase):

  def testConformerRMS(self):
    m1 = Chem.MultiConfMolFromSDF(os.path.join(RDConfig.RDCodeDir, "Chem/test_data/biaryl.sdf"),
                                  removeHs=False)

    m2 = Chem.MolFromSmiles('CNc(n2)nc(C)cc2Nc(cc34)ccc3[nH]nc4')
    m2.AddConformer(m1.GetConformer(id=1))

    # test that the prealigned flag is working
    rms1 = AllChem.GetConformerRMS(m1, 0, 1, prealigned=True)
    rms2 = AllChem.GetConformerRMS(m1, 0, 1, prealigned=False)
    self.assertTrue((rms1 > rms2))

    # test that RMS is the same as calculated by AlignMol()
    self.assertAlmostEqual(rms2, AllChem.GetBestRMS(m2, m1, 1, 0), 3)

    # the RMS with itself must be zero
    rms2 = AllChem.GetConformerRMS(m1, 0, 0, prealigned=True)
    self.assertAlmostEqual(rms2, 0.0, 4)

  def testConformerRMSMatrixAlignmentSubset(self):
    mol = Chem.MolFromSmiles('CCCC')
    # The first three atoms superimpose. The fourth rotates around the central bond.
    for offset, end in enumerate(((2., 1.4, 0.), (2., -1.4, 0.), (2., 0., 1.4))):
      conf = Chem.Conformer(4)
      for idx, (x, y, z) in enumerate(((-0.5, 1.4, 0.), (0., 0., 0.), (1.5, 0., 0.), end)):
        conf.SetAtomPosition(idx, (x + offset, y, z))
      mol.AddConformer(conf, assignId=True)

    atoms = [0, 1, 2]
    # Without alignment the selected atoms differ only by the translations.
    self.assertEqual(AllChem.GetConformerRMSMatrix(mol, atoms, prealigned=True), [1., 2., 1.])
    aligned = AllChem.GetConformerRMSMatrix(mol, atoms)
    for result in (aligned, AllChem.GetConformerRMSMatrix(mol, atoms, prealigned=True)):
      self.assertEqual(len(result), 3)
      for rms in result:
        self.assertAlmostEqual(rms, 0.)

    # Measuring all atoms on the same coordinates includes the moving fourth atom.
    expected = [1.4, 1.4 / 2**0.5, 1.4 / 2**0.5]
    for atomIds in (None, []):
      result = AllChem.GetConformerRMSMatrix(mol, atomIds, prealigned=True)
      self.assertEqual(len(result), len(expected))
      for rms, ref in zip(result, expected):
        self.assertAlmostEqual(rms, ref)

  def testConformerRMSMatrixScaledSubset(self):
    mol = Chem.MolFromSmiles('CCCC')
    for s in (1., 2., 3.):
      conf = Chem.Conformer(4)
      for idx, xyz in enumerate(((s, 0., 0.), (0., 0., 1.), (0., s, 0.), (-s, -s, 0.))):
        conf.SetAtomPosition(idx, xyz)
      mol.AddConformer(conf, assignId=True)

    # The selected triangle expands around its centroid. Between adjacent scales
    # the squared displacements are 1, 1 and 2, so its RMS is sqrt(4/3).
    rms = (4. / 3.)**0.5
    for prealigned in (False, True):
      result = AllChem.GetConformerRMSMatrix(mol, [0, 2, 3], prealigned)
      self.assertEqual(len(result), 3)
      for actual, expected in zip(result, (rms, 2. * rms, rms)):
        self.assertAlmostEqual(actual, expected)

  def testConformerRMSMatrix(self):
    m1 = Chem.MultiConfMolFromSDF(os.path.join(RDConfig.RDCodeDir, "Chem/test_data/biaryl.sdf"),
                                  removeHs=False)

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
    refFile = os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data', '1DWD_ligand.pdb')
    ref = Chem.MolFromSmiles(
      'NC(=[NH2+])c1ccc(C[C@@H](NC(=O)CNS(=O)(=O)c2ccc3ccccc3c2)C(=O)N2CCCCC2)cc1')
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
    self.assertEqual(len(weights), len(tors_list + tors_list_rings))
    weights = TorsionFingerprints.CalculateTorsionWeights(mol, 15, 14)
    self.assertAlmostEqual(weights[3], 1.0)
    self.assertRaises(ValueError, TorsionFingerprints.CalculateTorsionWeights, mol, 15, 3)

    # the torsion angles
    tors_list, tors_list_rings = TorsionFingerprints.CalculateTorsionLists(mol)
    torsions = TorsionFingerprints.CalculateTorsionAngles(mol, tors_list, tors_list_rings)
    self.assertEqual(len(weights), len(torsions))
    self.assertAlmostEqual(torsions[2][0][0], 232.5346, 4)

    # the torsion fingerprint deviation
    tfd = TorsionFingerprints.CalculateTFD(torsions, torsions)
    self.assertAlmostEqual(tfd, 0.0)
    refFile = os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data', '1PPC_ligand.pdb')
    mol2 = Chem.MolFromPDBFile(refFile)
    mol2 = AllChem.AssignBondOrdersFromTemplate(ref, mol2)
    torsions2 = TorsionFingerprints.CalculateTorsionAngles(mol2, tors_list, tors_list_rings)
    weights = TorsionFingerprints.CalculateTorsionWeights(mol)
    tfd = TorsionFingerprints.CalculateTFD(torsions, torsions2, weights=weights)
    self.assertAlmostEqual(tfd, 0.0691, 4)
    tfd = TorsionFingerprints.CalculateTFD(torsions, torsions2)
    self.assertAlmostEqual(tfd, 0.1115, 4)

    # the wrapper functions
    tfd = TorsionFingerprints.GetTFDBetweenMolecules(mol, mol2)
    self.assertAlmostEqual(tfd, 0.0691, 4)

    # exactly equivalent to the above since the mols each only have one conformer
    tfd = TorsionFingerprints.GetBestTFDBetweenMolecules(mol2, mol)
    self.assertAlmostEqual(tfd, 0.0691, 4)

    mol.AddConformer(mol2.GetConformer(), assignId=True)
    mol.AddConformer(mol2.GetConformer(), assignId=True)
    tfd = TorsionFingerprints.GetTFDBetweenConformers(mol, confIds1=[0], confIds2=[1, 2])
    self.assertEqual(len(tfd), 2)
    self.assertAlmostEqual(tfd[0], 0.0691, 4)

    tfdmat = TorsionFingerprints.GetTFDMatrix(mol)
    self.assertEqual(len(tfdmat), 3)

    tfd = TorsionFingerprints.GetBestTFDBetweenMolecules(mol2, mol)
    self.assertAlmostEqual(tfd, 0, 4)

  def testTorsionFingerprintsAtomReordering(self):
    # we use the xray structure from the paper (JCIM, 52, 1499, 2012): 1DWD
    refFile = os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data', '1DWD_ligand.pdb')
    ref = Chem.MolFromSmiles(
      'NC(=[NH2+])c1ccc(C[C@@H](NC(=O)CNS(=O)(=O)c2ccc3ccccc3c2)C(=O)N2CCCCC2)cc1')
    mol1 = Chem.MolFromPDBFile(refFile)
    mol1 = AllChem.AssignBondOrdersFromTemplate(ref, mol1)

    refFile = os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data', '1DWD_ligand_reordered.pdb')
    mol2 = Chem.MolFromPDBFile(refFile)
    mol2 = AllChem.AssignBondOrdersFromTemplate(ref, mol2)

    tfd = TorsionFingerprints.GetTFDBetweenMolecules(mol1, mol2)
    self.assertEqual(tfd, 0.0)

  def testTorsionFingerprintsColinearBonds(self):
    # test that single bonds adjacent to triple bonds are ignored
    mol = Chem.MolFromSmiles('CCC#CCC')
    tors_list, tors_list_rings = TorsionFingerprints.CalculateTorsionLists(
      mol, ignoreColinearBonds=True)
    self.assertEqual(len(tors_list), 0)
    weights = TorsionFingerprints.CalculateTorsionWeights(mol, ignoreColinearBonds=True)
    self.assertEqual(len(weights), 0)

    # test that they are not ignored, but alternative atoms searched for
    tors_list, tors_list_rings = TorsionFingerprints.CalculateTorsionLists(
      mol, ignoreColinearBonds=False)
    self.assertEqual(len(tors_list), 1)
    self.assertEqual(tors_list[0][0][0], (0, 1, 4, 5))
    weights = TorsionFingerprints.CalculateTorsionWeights(mol, ignoreColinearBonds=False)
    self.assertEqual(len(weights), 1)

    # test that single bonds adjacent to terminal triple bonds are always ignored
    mol = Chem.MolFromSmiles('C#CCC')
    tors_list, tors_list_rings = TorsionFingerprints.CalculateTorsionLists(
      mol, ignoreColinearBonds=True)
    self.assertEqual(len(tors_list), 0)
    tors_list, tors_list_rings = TorsionFingerprints.CalculateTorsionLists(
      mol, ignoreColinearBonds=False)
    self.assertEqual(len(tors_list), 0)

  def testGithub4720(self):
    # exceptions with highly-coordinated atoms
    mol = Chem.MolFromSmiles('S(F)(F)(F)(F)(Cl)c1ccccc1')
    tors_list, tors_list_rings = TorsionFingerprints.CalculateTorsionLists(mol)
    self.assertEqual(len(tors_list), 1)

  def testTorsionAngleLeargerThan14(self):
    # incorrect value from more than 15-membered ring
    mol = Chem.MolFromSmiles('C1' + 'C' * 13 + 'C1')
    tors_list, tors_list_rings = TorsionFingerprints.CalculateTorsionLists(mol)
    self.assertAlmostEqual(tors_list_rings[-1][1], 180.0, 4)

  def assertBondStereoRoundTrips(self, fname):
    path = os.path.join(RDConfig.RDCodeDir, 'Chem', 'test_data', fname)
    mol = Chem.MolFromMolFile(path)
    refSmiles = mol.GetProp("_Name")
    self.assertTrue(len(refSmiles) > 0)
    self.assertEqual(Chem.MolToSmiles(mol, isomericSmiles=True), refSmiles)

    # now test Chem.DetectBondStereoChemistry more directly by constructing the molecule from scratch
    oldconf = mol.GetConformer(0)
    newconf = Chem.Conformer(mol.GetNumAtoms())
    newmol = Chem.RWMol()

    for atm in mol.GetAtoms():
      ratm = Chem.Atom(atm.GetAtomicNum())
      ratm.SetFormalCharge(atm.GetFormalCharge())
      newmol.AddAtom(ratm)

      atomidx = atm.GetIdx()
      pos = oldconf.GetAtomPosition(atomidx)
      newconf.SetAtomPosition(atomidx, pos)

    for bnd in mol.GetBonds():
      newmol.AddBond(bnd.GetBeginAtomIdx(), bnd.GetEndAtomIdx(), Chem.BondType(bnd.GetBondType()))
    newmol.AddConformer(newconf)

    Chem.SanitizeMol(newmol)
    Chem.DetectBondStereoChemistry(newmol, newmol.GetConformer())

    # these aren't necessary for this specific test case, but are for
    # a more general conversion routine, so would like to see them
    # tested eventually
    # Chem.AssignAtomChiralTagsFromStructure(newmol)
    # Chem.AssignStereochemistry(newmol)

    self.assertEqual(Chem.MolToSmiles(newmol, isomericSmiles=True), refSmiles)

  def testDetectBondStereoChemistry(self):
    self.assertBondStereoRoundTrips('cis.sdf')
    self.assertBondStereoRoundTrips('trans.sdf')


if __name__ == '__main__':
  unittest.main()
