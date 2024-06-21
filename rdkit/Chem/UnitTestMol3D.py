"""unit testing code for 3D stuff

"""
import os
import random
import sys
import unittest

from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, TorsionFingerprints


class TestCase(unittest.TestCase):

  def testConformerRMS(self):
    m1 = Chem.MolFromSmiles('CNc(n2)nc(C)cc2Nc(cc34)ccc3[nH]nc4')
    cids = AllChem.EmbedMultipleConfs(m1, 2)

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

  def testConformerRMSMatrix(self):
    m1 = Chem.MolFromSmiles('CNc(n2)nc(C)cc2Nc(cc34)ccc3[nH]nc4')
    cids = AllChem.EmbedMultipleConfs(m1, 3)

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

  def testEnumerateStereoisomersBasic(self):
    mol = Chem.MolFromSmiles('CC(F)=CC(Cl)C')
    smiles = set(
      Chem.MolToSmiles(i, isomericSmiles=True) for i in AllChem.EnumerateStereoisomers(mol))
    self.assertEqual(len(smiles), 4)

  def testEnumerateStereoisomersLargeRandomSample(self):
    # near max number of stereo centers allowed
    mol = Chem.MolFromSmiles('CC(F)=CC(Cl)C' * 31)
    smiles = set(
      Chem.MolToSmiles(i, isomericSmiles=True) for i in AllChem.EnumerateStereoisomers(mol))
    self.assertEqual(len(smiles), 1024)

  def testEnumerateStereoisomersWithCrazyNumberOfCenters(self):
    # insanely large numbers of isomers aren't a problem
    mol = Chem.MolFromSmiles('CC(F)=CC(Cl)C' * 101)
    opts = AllChem.StereoEnumerationOptions(rand=None, maxIsomers=13)
    smiles = set(
      Chem.MolToSmiles(i, isomericSmiles=True) for i in AllChem.EnumerateStereoisomers(mol, opts))
    self.assertEqual(len(smiles), 13)

  def testEnumerateStereoisomersRandomSamplingShouldBeDeterministicAndPortable(self):
    mol = Chem.MolFromSmiles('CC(F)=CC(Cl)C=C(Br)C(I)N')
    opts = AllChem.StereoEnumerationOptions(maxIsomers=2)
    smiles = set(
      Chem.MolToSmiles(i, isomericSmiles=True) for i in AllChem.EnumerateStereoisomers(mol, opts))
    expected = set(
      ['C/C(F)=C\\[C@H](Cl)/C=C(/Br)[C@H](N)I', 'C/C(F)=C/[C@@H](Cl)/C=C(/Br)[C@H](N)I'])
    self.assertEqual(smiles, expected)

  def testEnumerateStereoisomersMaxIsomersShouldBeReturnedEvenWithTryEmbedding(self):
    m = Chem.MolFromSmiles('BrC=CC1OC(C2)(F)C2(Cl)C1')
    opts = AllChem.StereoEnumerationOptions(tryEmbedding=True, maxIsomers=8)
    isomers = set()
    for x in AllChem.EnumerateStereoisomers(m, options=opts):
      isomers.add(Chem.MolToSmiles(x, isomericSmiles=True))
    self.assertEqual(len(isomers), 8)

  def testEnumerateStereoisomersTryEmbeddingShouldNotInfiniteLoopWhenMaxIsomersIsLargerThanActual(
      self):
    m = Chem.MolFromSmiles('BrC=CC1OC(C2)(F)C2(Cl)C1')
    opts = AllChem.StereoEnumerationOptions(tryEmbedding=True, maxIsomers=1024)
    isomers = set()
    for x in AllChem.EnumerateStereoisomers(m, options=opts):
      isomers.add(Chem.MolToSmiles(x, isomericSmiles=True))
    self.assertEqual(len(isomers), 8)

  def testEnumerateStereoisomersRandomSeeding(self):
    opts = AllChem.StereoEnumerationOptions(rand=None, maxIsomers=3)
    mol = Chem.MolFromSmiles('CC(F)=CC(Cl)C')
    smiles = set(
      Chem.MolToSmiles(i, isomericSmiles=True) for i in AllChem.EnumerateStereoisomers(mol, opts))
    self.assertEqual(smiles,
                     set(['C/C(F)=C/[C@@H](C)Cl', 'C/C(F)=C\\[C@H](C)Cl', 'C/C(F)=C\\[C@@H](C)Cl']))

    opts = AllChem.StereoEnumerationOptions(rand=0xDEADBEEF)
    mol = Chem.MolFromSmiles('c1ccc2c(c1)C(=O)N(C2=O)C3CCC(=O)NC3=O')
    smiles = set(
      Chem.MolToSmiles(i, isomericSmiles=True) for i in AllChem.EnumerateStereoisomers(mol, opts))
    self.assertEqual(
      smiles,
      set(['O=C1CC[C@@H](N2C(=O)c3ccccc3C2=O)C(=O)N1', 'O=C1CC[C@H](N2C(=O)c3ccccc3C2=O)C(=O)N1']))

    class DeterministicRandom(random.Random):

      def __init__(self):
        random.Random.__init__(self)
        self.count = 0

      def getrandbits(self, n_bits):
        c = self.count
        self.count += 1
        return c

    # start with the non-seeded version, just to make sure we get everything
    mol = Chem.MolFromSmiles('CCCC(=C(CCl)C(C)CBr)[C@H](F)C(C)C')
    opts = AllChem.StereoEnumerationOptions()
    smiles = list(
      Chem.MolToSmiles(i, isomericSmiles=True) for i in AllChem.EnumerateStereoisomers(mol))
    self.assertEqual(smiles, [
      'CCC/C(=C(/CCl)[C@H](C)CBr)[C@H](F)C(C)C', 'CCC/C(=C(/CCl)[C@@H](C)CBr)[C@H](F)C(C)C',
      'CCC/C(=C(\\CCl)[C@H](C)CBr)[C@H](F)C(C)C', 'CCC/C(=C(\\CCl)[C@@H](C)CBr)[C@H](F)C(C)C'
    ])
    rand = DeterministicRandom()
    opts = AllChem.StereoEnumerationOptions(rand=rand, maxIsomers=3)
    smiles = list(
      sorted(
        Chem.MolToSmiles(i, isomericSmiles=True)
        for i in AllChem.EnumerateStereoisomers(mol, opts)))
    self.assertEqual(smiles, [
      'CCC/C(=C(/CCl)[C@@H](C)CBr)[C@H](F)C(C)C', 'CCC/C(=C(/CCl)[C@H](C)CBr)[C@H](F)C(C)C',
      'CCC/C(=C(\\CCl)[C@H](C)CBr)[C@H](F)C(C)C'
    ])

  def testEnumerateStereoisomersOnlyUnassigned(self):
    # shouldn't enumerate anything
    fully_assigned = Chem.MolFromSmiles('C/C(F)=C/[C@@H](C)Cl')
    smiles = set(
      Chem.MolToSmiles(i, isomericSmiles=True)
      for i in AllChem.EnumerateStereoisomers(fully_assigned))
    self.assertEqual(smiles, set(['C/C(F)=C/[C@@H](C)Cl']))

    # should only enumerate the bond stereo
    partially_assigned = Chem.MolFromSmiles('CC(F)=C[C@@H](C)Cl')
    smiles = set(
      Chem.MolToSmiles(i, isomericSmiles=True)
      for i in AllChem.EnumerateStereoisomers(partially_assigned))
    self.assertEqual(smiles, set(['C/C(F)=C/[C@@H](C)Cl', 'C/C(F)=C\\[C@@H](C)Cl']))

    # should enumerate everything
    opts = AllChem.StereoEnumerationOptions(onlyUnassigned=False)
    smiles = set(
      Chem.MolToSmiles(i, isomericSmiles=True)
      for i in AllChem.EnumerateStereoisomers(fully_assigned, opts))
    self.assertEqual(
      smiles,
      set([
        'C/C(F)=C\\[C@@H](C)Cl',
        'C/C(F)=C\\[C@H](C)Cl',
        'C/C(F)=C/[C@H](C)Cl',
        'C/C(F)=C/[C@@H](C)Cl',
      ]))

  def testEnumerateStereoisomersOnlyUnique(self):
    mol = Chem.MolFromSmiles('FC(Cl)C(Cl)F')
    opts = AllChem.StereoEnumerationOptions(unique=False)
    smiles = [
      Chem.MolToSmiles(i, isomericSmiles=True) for i in AllChem.EnumerateStereoisomers(mol, opts)
    ]
    self.assertEqual(len(smiles), 4)
    self.assertEqual(len(set(smiles)), 3)

    mol = Chem.MolFromSmiles('FC(Cl)C(Cl)F')
    opts = AllChem.StereoEnumerationOptions(unique=True)
    smiles = set(
      Chem.MolToSmiles(i, isomericSmiles=True) for i in AllChem.EnumerateStereoisomers(mol, opts))
    self.assertEqual(
      smiles, set(['F[C@@H](Cl)[C@@H](F)Cl', 'F[C@H](Cl)[C@@H](F)Cl', 'F[C@H](Cl)[C@H](F)Cl']))

    mol = Chem.MolFromSmiles('CC=CC=CC')
    opts = AllChem.StereoEnumerationOptions(unique=False)
    smiles = [
      Chem.MolToSmiles(i, isomericSmiles=True) for i in AllChem.EnumerateStereoisomers(mol, opts)
    ]
    self.assertEqual(len(smiles), 4)
    self.assertEqual(len(set(smiles)), 3)

    mol = Chem.MolFromSmiles('CC=CC=CC')
    opts = AllChem.StereoEnumerationOptions(unique=True)
    smiles = set(
      Chem.MolToSmiles(i, isomericSmiles=True) for i in AllChem.EnumerateStereoisomers(mol, opts))
    self.assertEqual(smiles, set(['C/C=C/C=C/C', 'C/C=C\\C=C\\C', 'C/C=C\\C=C/C']))

    mol = Chem.MolFromSmiles('FC(Cl)C=CC=CC(F)Cl')
    opts = AllChem.StereoEnumerationOptions(unique=False)
    smiles = [
      Chem.MolToSmiles(i, isomericSmiles=True) for i in AllChem.EnumerateStereoisomers(mol, opts)
    ]
    self.assertEqual(len(smiles), 16)
    self.assertEqual(len(set(smiles)), 10)

    mol = Chem.MolFromSmiles('FC(Cl)C=CC=CC(F)Cl')
    opts = AllChem.StereoEnumerationOptions(unique=True)
    smiles = set(
      sorted(
        Chem.MolToSmiles(i, isomericSmiles=True)
        for i in AllChem.EnumerateStereoisomers(mol, opts)))
    self.assertEqual(
      smiles,
      set(
        sorted([
          'F[C@H](Cl)/C=C/C=C/[C@@H](F)Cl', 'F[C@H](Cl)/C=C/C=C/[C@H](F)Cl',
          'F[C@@H](Cl)/C=C/C=C/[C@@H](F)Cl', 'F[C@@H](Cl)/C=C\\C=C\\[C@@H](F)Cl',
          'F[C@H](Cl)/C=C\\C=C/[C@@H](F)Cl', 'F[C@H](Cl)/C=C/C=C\\[C@@H](F)Cl',
          'F[C@H](Cl)/C=C\\C=C\\[C@H](F)Cl', 'F[C@@H](Cl)/C=C\\C=C/[C@@H](F)Cl',
          'F[C@H](Cl)/C=C\\C=C/[C@H](F)Cl', 'F[C@H](Cl)/C=C\\C=C\\[C@@H](F)Cl'
        ])))

  def testEnumerateStereoisomersOnlyEnhancedStereo(self):
    rdbase = os.environ["RDBASE"]
    filename = os.path.join(rdbase, 'Code/GraphMol/FileParsers/test_data/two_centers_or.mol')
    mol = Chem.MolFromMolFile(filename)
    smiles = set(Chem.MolToSmiles(m) for m in AllChem.EnumerateStereoisomers(mol))
    # switches the centers linked by an "OR", but not the absolute group
    self.assertEqual(smiles, {r'C[C@H]([C@@H](C)F)[C@@H](C)Br', r'C[C@@H]([C@H](C)F)[C@@H](C)Br'})

    # we need the SMILES without the influence of the stereo groups:
    mol = Chem.RWMol(mol)
    mol.SetStereoGroups([])
    original_smiles = Chem.MolToSmiles(mol)
    self.assertIn(original_smiles, smiles)

  def testNoExtrasIfEnumeratingAllWithEnhancedStereo(self):
    """
    If the onlyUnassigned option is False, make sure that enhanced stereo
    groups aren't double-counted.
    """
    rdbase = os.environ["RDBASE"]
    filename = os.path.join(rdbase, 'Code/GraphMol/FileParsers/test_data/two_centers_or.mol')
    mol = Chem.MolFromMolFile(filename)

    opts = AllChem.StereoEnumerationOptions(onlyUnassigned=False, unique=False)
    smiles = [Chem.MolToSmiles(m) for m in AllChem.EnumerateStereoisomers(mol, opts)]
    self.assertEqual(len(smiles), len(set(smiles)))
    self.assertEqual(len(smiles), 2**3)

  def testIssue2890(self):
    mol = Chem.MolFromSmiles('CC=CC')
    mol.GetBondWithIdx(1).SetStereo(Chem.rdchem.BondStereo.STEREOANY)

    self.assertEqual(len(list(AllChem.EnumerateStereoisomers(mol))), 2)

  def testIssue3231(self):
    mol = Chem.MolFromSmiles(
      'C[C@H](OC1=C(N)N=CC(C2=CN(C3C[C@H](C)NCC3)N=C2)=C1)C4=C(Cl)C=CC(F)=C4Cl')
    Chem.AssignStereochemistry(mol, force=True, flagPossibleStereoCenters=True)
    l = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    self.assertEqual(l, [(1, 'S'), (12, '?'), (14, 'S')])
    enumsi_opt = AllChem.StereoEnumerationOptions(maxIsomers=20, onlyUnassigned=False)
    isomers = list(AllChem.EnumerateStereoisomers(mol, enumsi_opt))
    chi_cents = []
    for iso in isomers:
      Chem.AssignStereochemistry(iso)
      chi_cents.append(Chem.FindMolChiralCenters(iso, includeUnassigned=True))
    self.assertEqual(sorted(chi_cents),
                     [[(1, 'R'), (12, 'R'),
                       (14, 'R')], [(1, 'R'), (12, 'R'),
                                    (14, 'S')], [(1, 'R'), (12, 'S'),
                                                 (14, 'R')], [(1, 'R'), (12, 'S'), (14, 'S')],
                      [(1, 'S'), (12, 'R'),
                       (14, 'R')], [(1, 'S'), (12, 'R'),
                                    (14, 'S')], [(1, 'S'), (12, 'S'),
                                                 (14, 'R')], [(1, 'S'), (12, 'S'), (14, 'S')]])

  def testIssue3505(self):
    m = Chem.MolFromSmiles('CCC(C)Br')
    mols = list(AllChem.EnumerateStereoisomers(m))
    self.assertEqual(len(mols), 2)
    for mol in mols:
      at = mol.GetAtomWithIdx(2)
      self.assertIn(at.GetChiralTag(),
                    [Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW])
      self.assertTrue(at.HasProp("_ChiralityPossible"))

  def testEnumerateEitherDoubleStereo(self):
    """ EnumerateStereoisomers from MOL with explicit either cis/trans bond """
    rdbase = os.environ["RDBASE"]
    filename = os.path.join(rdbase, 'Code/GraphMol/FileParsers/test_data/simple_either.mol')
    mol = Chem.MolFromMolFile(filename)
    smiles = [Chem.MolToSmiles(m) for m in AllChem.EnumerateStereoisomers(mol)]
    self.assertEqual(set(smiles), {"C/C=C/C", "C/C=C\\C"})

  def testTryEmbeddingManyChirals(self):
    smiles = "C1" + "C(Cl)(Br)" * 40 + "C1"
    mol = Chem.MolFromSmiles(smiles)
    opts = AllChem.StereoEnumerationOptions(tryEmbedding=True, maxIsomers=2)
    self.assertEqual(len(list(AllChem.EnumerateStereoisomers(mol, options=opts))), 2)

  def testGithub6045(self):
    mol = Chem.MolFromSmiles('O[C@H](Br)[C@H](F)C |&1:1,3|')
    prods = list(AllChem.EnumerateStereoisomers(mol))
    self.assertEqual(len(prods), 2)
    for prod in prods:
      self.assertEqual(len(prod.GetStereoGroups()), 0)

  def testGithub7516(self):
    m = Chem.MolFromSmiles('CC1CC(C)C1')
    sis = list(AllChem.EnumerateStereoisomers(m))
    self.assertEqual(len(sis), 2)
    self.assertEqual(Chem.MolToSmiles(sis[0]), 'C[C@H]1C[C@@H](C)C1')
    self.assertEqual(Chem.MolToSmiles(sis[1]), 'C[C@H]1C[C@H](C)C1')

    m = Chem.MolFromSmiles('COC(=O)C1CC(NC(N)=O)C1')
    sis = list(AllChem.EnumerateStereoisomers(m))
    self.assertEqual(len(sis), 2)

    m = Chem.MolFromSmiles('O=C(NC1CC2[NH+](C(C1)CC2)Cc3ccccc3)N')
    sis = list(AllChem.EnumerateStereoisomers(m))
    self.assertEqual(len(sis), 8)


if __name__ == '__main__':
  unittest.main()
