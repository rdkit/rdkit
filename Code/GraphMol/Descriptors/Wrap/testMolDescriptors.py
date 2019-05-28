# $Id$
#

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors as rdMD, Descriptors
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit import RDConfig
from rdkit.Geometry import rdGeometry as rdG
import unittest


def feq(v1, v2, tol=1.e-4):
  return abs(v1 - v2) < tol


class TestCase(unittest.TestCase):

  def setUp(self):
    pass

  def testAtomPairTypes(self):
    params = rdMD.AtomPairsParameters
    mol = Chem.MolFromSmiles("C=C")
    self.assertTrue(rdMD.GetAtomPairAtomCode(mol.GetAtomWithIdx(0))==\
                    rdMD.GetAtomPairAtomCode(mol.GetAtomWithIdx(1)))
    self.assertTrue(rdMD.GetAtomPairAtomCode(mol.GetAtomWithIdx(0))==\
                    1 | (1 | 1<<params.numPiBits)<<params.numBranchBits)

    mol = Chem.MolFromSmiles("C#CO")
    self.assertTrue(rdMD.GetAtomPairAtomCode(mol.GetAtomWithIdx(0))!=\
                    rdMD.GetAtomPairAtomCode(mol.GetAtomWithIdx(1)))
    self.assertTrue(rdMD.GetAtomPairAtomCode(mol.GetAtomWithIdx(0))==\
                    1 | (2 | 1<<params.numPiBits)<<params.numBranchBits)
    self.assertTrue(rdMD.GetAtomPairAtomCode(mol.GetAtomWithIdx(1))==\
                    2 | (2 | 1<<params.numPiBits)<<params.numBranchBits)
    self.assertTrue(rdMD.GetAtomPairAtomCode(mol.GetAtomWithIdx(2))==\
                    1 | (0 | 3<<params.numPiBits)<<params.numBranchBits)
    self.assertTrue(rdMD.GetAtomPairAtomCode(mol.GetAtomWithIdx(1),1)==\
                    1 | (2 | 1<<params.numPiBits)<<params.numBranchBits)
    self.assertTrue(rdMD.GetAtomPairAtomCode(mol.GetAtomWithIdx(1),2)==\
                    0 | (2 | 1<<params.numPiBits)<<params.numBranchBits)

  def testAtomPairTypesChirality(self):
    mols = [Chem.MolFromSmiles(x) for x in ("CC(F)Cl", "C[C@@H](F)Cl", "C[C@H](F)Cl")]
    self.assertEqual(
      rdMD.GetAtomPairAtomCode(mols[0].GetAtomWithIdx(1)),
      rdMD.GetAtomPairAtomCode(mols[1].GetAtomWithIdx(1)))
    self.assertEqual(
      rdMD.GetAtomPairAtomCode(mols[0].GetAtomWithIdx(1)),
      rdMD.GetAtomPairAtomCode(mols[2].GetAtomWithIdx(1)))
    self.assertEqual(
      rdMD.GetAtomPairAtomCode(mols[0].GetAtomWithIdx(1), includeChirality=True),
      rdMD.GetAtomPairAtomCode(mols[0].GetAtomWithIdx(1)))
    self.assertNotEqual(
      rdMD.GetAtomPairAtomCode(mols[0].GetAtomWithIdx(1), includeChirality=True),
      rdMD.GetAtomPairAtomCode(mols[1].GetAtomWithIdx(1), includeChirality=True))
    self.assertNotEqual(
      rdMD.GetAtomPairAtomCode(mols[0].GetAtomWithIdx(1), includeChirality=True),
      rdMD.GetAtomPairAtomCode(mols[2].GetAtomWithIdx(1), includeChirality=True))
    self.assertNotEqual(
      rdMD.GetAtomPairAtomCode(mols[1].GetAtomWithIdx(1), includeChirality=True),
      rdMD.GetAtomPairAtomCode(mols[2].GetAtomWithIdx(1), includeChirality=True))

    fps = [rdMD.GetAtomPairFingerprint(x) for x in mols]
    chiralFps = [rdMD.GetAtomPairFingerprint(x, includeChirality=True) for x in mols]
    for mol, fp, cfp in zip(mols, fps, chiralFps):
      ac0 = rdMD.GetAtomPairAtomCode(mol.GetAtomWithIdx(0))
      ac1 = rdMD.GetAtomPairAtomCode(mol.GetAtomWithIdx(1))
      self.assertTrue(rdMD.GetAtomPairCode(ac0, ac1, 1) in fp.GetNonzeroElements())
      ac0 = rdMD.GetAtomPairAtomCode(mol.GetAtomWithIdx(0), includeChirality=True)
      ac1 = rdMD.GetAtomPairAtomCode(mol.GetAtomWithIdx(1), includeChirality=True)
      self.assertFalse(
        rdMD.GetAtomPairCode(ac0, ac1, 1, includeChirality=True) in fp.GetNonzeroElements())
      self.assertTrue(
        rdMD.GetAtomPairCode(ac0, ac1, 1, includeChirality=True) in cfp.GetNonzeroElements())

  def testAtomPairs(self):
    m = Chem.MolFromSmiles('CCC')
    fp1 = rdMD.GetAtomPairFingerprint(m)
    fp2 = rdMD.GetAtomPairFingerprint(m, minLength=1, maxLength=2)
    nz1 = fp1.GetNonzeroElements()
    self.assertEqual(len(nz1), 2)
    nz2 = fp2.GetNonzeroElements()
    self.assertEqual(len(nz2), 2)

    fp2 = rdMD.GetAtomPairFingerprint(m, minLength=1, maxLength=1)
    nz2 = fp2.GetNonzeroElements()
    self.assertEqual(len(nz2), 1)

  def testHashedAtomPairs(self):
    m = Chem.MolFromSmiles('c1ccccc1')
    fp1 = rdMD.GetHashedAtomPairFingerprint(m, 2048)
    fp2 = rdMD.GetHashedAtomPairFingerprint(m, 2048, 1, 3)
    self.assertTrue(fp1 == fp2)
    fp2 = rdMD.GetHashedAtomPairFingerprint(m, 2048, 1, 2)
    sim = DataStructs.DiceSimilarity(fp1, fp2)
    self.assertTrue(sim > 0.0 and sim < 1.0)

    m = Chem.MolFromSmiles('c1ccccn1')
    fp2 = rdMD.GetHashedAtomPairFingerprint(m, 2048)
    sim = DataStructs.DiceSimilarity(fp1, fp2)
    self.assertTrue(sim > 0.0 and sim < 1.0)

    m = Chem.MolFromSmiles('c1ccccc1')
    fp1 = rdMD.GetHashedAtomPairFingerprintAsBitVect(m, 2048)
    m = Chem.MolFromSmiles('c1ccccn1')
    fp2 = rdMD.GetHashedAtomPairFingerprintAsBitVect(m, 2048)
    sim = DataStructs.DiceSimilarity(fp1, fp2)
    self.assertTrue(sim > 0.0 and sim < 1.0)

  def testRootedAtomPairs(self):
    m = Chem.MolFromSmiles('Oc1ccccc1')
    fp1 = rdMD.GetAtomPairFingerprint(m)
    fp2 = rdMD.GetAtomPairFingerprint(m, fromAtoms=(0, ))
    nz1 = fp1.GetNonzeroElements()
    nz2 = fp2.GetNonzeroElements()
    for k, v in nz2.items():
      self.assertTrue(v <= nz1[k])

  def testTopologicalTorsions(self):
    mol = Chem.MolFromSmiles("CC")
    fp = rdMD.GetTopologicalTorsionFingerprint(mol)
    self.assertTrue(fp.GetTotalVal() == 0)

    mol = Chem.MolFromSmiles("CCCC")
    fp = rdMD.GetTopologicalTorsionFingerprint(mol)
    self.assertTrue(fp.GetTotalVal() == 1)
    fp = rdMD.GetTopologicalTorsionFingerprint(mol, 3)
    self.assertTrue(fp.GetTotalVal() == 2)

    mol = Chem.MolFromSmiles("CCCO")
    fp = rdMD.GetTopologicalTorsionFingerprint(mol)
    self.assertTrue(fp.GetTotalVal() == 1)
    fp = rdMD.GetTopologicalTorsionFingerprint(mol, 3)
    self.assertTrue(fp.GetTotalVal() == 2)

    mol = Chem.MolFromSmiles("CCCCCCCCCCC")
    fp = rdMD.GetTopologicalTorsionFingerprint(mol, 7)
    self.assertRaises(ValueError, lambda: rdMD.GetTopologicalTorsionFingerprint(mol, 8))

  def testHashedTopologicalTorsions(self):
    mol = Chem.MolFromSmiles("c1ncccc1")
    fp1 = rdMD.GetHashedTopologicalTorsionFingerprint(mol)
    mol = Chem.MolFromSmiles("n1ccccc1")
    fp2 = rdMD.GetHashedTopologicalTorsionFingerprint(mol)
    self.assertEqual(DataStructs.DiceSimilarity(fp1, fp2), 1.0)

  def testRootedTorsions(self):
    m = Chem.MolFromSmiles('Oc1ccccc1')
    fp1 = rdMD.GetTopologicalTorsionFingerprint(m)
    fp2 = rdMD.GetTopologicalTorsionFingerprint(m, fromAtoms=(0, ))
    nz1 = fp1.GetNonzeroElements()
    nz2 = fp2.GetNonzeroElements()
    for k, v in nz2.items():
      self.assertTrue(v <= nz1[k])

    m = Chem.MolFromSmiles('COCC')
    fp1 = rdMD.GetTopologicalTorsionFingerprint(m)
    self.assertEqual(len(fp1.GetNonzeroElements()), 1)
    fp1 = rdMD.GetTopologicalTorsionFingerprint(m, fromAtoms=(0, ))
    self.assertEqual(len(fp1.GetNonzeroElements()), 1)
    fp1 = rdMD.GetTopologicalTorsionFingerprint(m, fromAtoms=(1, ))
    self.assertEqual(len(fp1.GetNonzeroElements()), 0)

  def testMorganFingerprints(self):
    mol = Chem.MolFromSmiles('CC(F)(Cl)C(F)(Cl)C')
    fp = rdMD.GetMorganFingerprint(mol, 0)
    self.assertTrue(len(fp.GetNonzeroElements()) == 4)

    mol = Chem.MolFromSmiles('CC')
    fp = rdMD.GetMorganFingerprint(mol, 0)
    self.assertTrue(len(fp.GetNonzeroElements()) == 1)
    self.assertTrue(list(fp.GetNonzeroElements().values())[0] == 2)
    fp = rdMD.GetMorganFingerprint(mol, 0, useCounts=False)
    self.assertTrue(len(fp.GetNonzeroElements()) == 1)
    self.assertTrue(list(fp.GetNonzeroElements().values())[0] == 1)

    mol = Chem.MolFromSmiles('CC(F)(Cl)C(F)(Cl)C')
    fp = rdMD.GetHashedMorganFingerprint(mol, 0)
    self.assertTrue(len(fp.GetNonzeroElements()) == 4)
    fp = rdMD.GetMorganFingerprint(mol, 1)
    self.assertTrue(len(fp.GetNonzeroElements()) == 8)
    fp = rdMD.GetHashedMorganFingerprint(mol, 1)
    self.assertTrue(len(fp.GetNonzeroElements()) == 8)
    fp = rdMD.GetMorganFingerprint(mol, 2)
    self.assertTrue(len(fp.GetNonzeroElements()) == 9)

    mol = Chem.MolFromSmiles('CC(F)(Cl)[C@](F)(Cl)C')
    fp = rdMD.GetMorganFingerprint(mol, 0)
    self.assertTrue(len(fp.GetNonzeroElements()) == 4)
    fp = rdMD.GetMorganFingerprint(mol, 1)
    self.assertTrue(len(fp.GetNonzeroElements()) == 8)
    fp = rdMD.GetMorganFingerprint(mol, 2)
    self.assertTrue(len(fp.GetNonzeroElements()) == 9)
    fp = rdMD.GetMorganFingerprint(mol, 0, useChirality=True)
    self.assertTrue(len(fp.GetNonzeroElements()) == 4)
    fp = rdMD.GetMorganFingerprint(mol, 1, useChirality=True)
    self.assertTrue(len(fp.GetNonzeroElements()) == 9)
    fp = rdMD.GetMorganFingerprint(mol, 2, useChirality=True)
    self.assertTrue(len(fp.GetNonzeroElements()) == 10)

    mol = Chem.MolFromSmiles('CCCCC')
    fp = rdMD.GetMorganFingerprint(mol, 0, fromAtoms=(0, ))
    self.assertTrue(len(fp.GetNonzeroElements()) == 1)

    mol = Chem.MolFromSmiles('CC1CC1')
    vs1 = rdMD.GetConnectivityInvariants(mol)
    self.assertEqual(len(vs1), mol.GetNumAtoms())
    fp1 = rdMD.GetMorganFingerprint(mol, 2, invariants=vs1)
    fp2 = rdMD.GetMorganFingerprint(mol, 2)
    self.assertEqual(fp1, fp2)

    vs2 = rdMD.GetConnectivityInvariants(mol, False)
    self.assertEqual(len(vs2), mol.GetNumAtoms())
    self.assertNotEqual(vs1, vs2)
    fp1 = rdMD.GetMorganFingerprint(mol, 2, invariants=vs2)
    self.assertNotEqual(fp1, fp2)

    mol = Chem.MolFromSmiles('Cc1ccccc1')
    vs1 = rdMD.GetFeatureInvariants(mol)
    self.assertEqual(len(vs1), mol.GetNumAtoms())
    self.assertEqual(vs1[0], 0)
    self.assertNotEqual(vs1[1], 0)
    self.assertEqual(vs1[1], vs1[2])
    self.assertEqual(vs1[1], vs1[3])
    self.assertEqual(vs1[1], vs1[4])

    mol = Chem.MolFromSmiles('FCCCl')
    vs1 = rdMD.GetFeatureInvariants(mol)
    self.assertEqual(len(vs1), mol.GetNumAtoms())
    self.assertEqual(vs1[1], 0)
    self.assertEqual(vs1[2], 0)
    self.assertNotEqual(vs1[0], 0)
    self.assertEqual(vs1[0], vs1[3])

    fp1 = rdMD.GetMorganFingerprint(mol, 0, invariants=vs1)
    fp2 = rdMD.GetMorganFingerprint(mol, 0, useFeatures=True)
    self.assertEqual(fp1, fp2)

  def testCrippen(self):
    mol = Chem.MolFromSmiles("n1ccccc1CO")
    contribs = rdMD._CalcCrippenContribs(mol)
    self.assertEqual(len(contribs), mol.GetNumAtoms())

    ts = [0] * mol.GetNumAtoms()
    contribs = rdMD._CalcCrippenContribs(mol, force=True, atomTypes=ts)
    self.assertEqual(ts, [59, 25, 25, 25, 25, 28, 17, 69])

    ls = [''] * mol.GetNumAtoms()
    contribs = rdMD._CalcCrippenContribs(mol, force=True, atomTypeLabels=ls)
    self.assertEqual(ls, ['N11', 'C18', 'C18', 'C18', 'C18', 'C21', 'C10', 'O2'])

  def testUSR(self):
    mol = Chem.MolFromSmiles("CC")
    AllChem.Compute2DCoords(mol)
    self.assertRaises(ValueError, lambda: rdMD.GetUSR(mol))
    mol = Chem.MolFromSmiles("C1CCCCC1")
    mol = Chem.AddHs(mol)
    self.assertRaises(ValueError, lambda: rdMD.GetUSR(mol))
    AllChem.Compute2DCoords(mol)
    usr = rdMD.GetUSR(mol)
    self.assertEqual(len(usr), 12)

    self.assertRaises(ValueError, lambda: rdMD.GetUSRDistributions([]))

    conf = mol.GetConformer()
    coords = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
    dist = rdMD.GetUSRDistributions(coords)
    self.assertEqual(len(dist), 4)
    self.assertEqual(len(dist[0]), mol.GetNumAtoms())
    self.assertRaises(ValueError, lambda: rdMD.GetUSRFromDistributions([]))
    usr2 = rdMD.GetUSRFromDistributions(dist)
    self.assertEqual(usr, usr2)

    self.assertRaises(ValueError, lambda: rdMD.GetUSRDistributionsFromPoints(coords, []))
    p = []
    dist = rdMD.GetUSRDistributions(coords, p)
    self.assertRaises(ValueError, lambda: rdMD.GetUSRDistributionsFromPoints([], p))
    dist2 = rdMD.GetUSRDistributionsFromPoints(coords, p)
    usr2 = rdMD.GetUSRFromDistributions(dist2)
    self.assertEqual(usr, usr2)

    mol2 = Chem.MolFromSmiles("C1CCCCC1")
    mol2 = Chem.AddHs(mol2)
    AllChem.Compute2DCoords(mol2)
    usr2 = rdMD.GetUSR(mol2)
    self.assertRaises(ValueError, lambda: rdMD.GetUSRScore(usr, usr2[:2]))
    self.assertEqual(rdMD.GetUSRScore(usr, usr2), 1.0)

    m1 = [4.44, 2.98, 1.04, 4.55, 4.70, 0.23, 8.30, 16.69, -22.97, 7.37, 15.64, 0.51]
    m2 = [4.39, 3.11, 1.36, 4.50, 4.44, 0.09, 8.34, 16.78, -23.20, 7.15, 16.52, 0.13]
    self.assertAlmostEqual(rdMD.GetUSRScore(m1, m2), 0.812, 2)

  def testUSRCAT(self):
    mol = Chem.MolFromSmiles("CC")
    AllChem.Compute2DCoords(mol)
    self.assertRaises(ValueError, lambda: rdMD.GetUSRCAT(mol))
    mol = Chem.MolFromSmiles("C1CCCCC1")
    mol = Chem.AddHs(mol)
    self.assertRaises(ValueError, lambda: rdMD.GetUSRCAT(mol))
    AllChem.Compute2DCoords(mol)
    usr = rdMD.GetUSRCAT(mol)
    self.assertEqual(len(usr), 60)
    self.assertRaises(ValueError, lambda: rdMD.GetUSRCAT(mol, atomSelections=[]))
    atoms = [[1, 2, 3, 4, 5, 6], []]
    usr2 = rdMD.GetUSRCAT(mol, atomSelections=atoms)
    self.assertEqual(len(usr2), 36)
    atoms = [[1, 2, 3, 4, 5, 6], [], [], []]
    usr2 = rdMD.GetUSRCAT(mol, atomSelections=atoms)
    self.assertEqual(len(usr2), 60)
    self.assertEqual(rdMD.GetUSRScore(usr, usr2, weights=[1.0, 1.0, 1.0, 1.0, 1.0]), 1.0)

  def testMolWt(self):
    mol = Chem.MolFromSmiles("C")
    amw = rdMD._CalcMolWt(mol)
    self.assertTrue(feq(amw, 16.043, .001))
    amw = rdMD._CalcMolWt(mol, True)
    self.assertTrue(feq(amw, 12.011, .001))
    mol2 = Chem.AddHs(mol)
    amw = rdMD._CalcMolWt(mol2)
    self.assertTrue(feq(amw, 16.043, .001))
    amw = rdMD._CalcMolWt(mol2, True)
    self.assertTrue(feq(amw, 12.011, .001))

    mol = Chem.MolFromSmiles("C")
    amw = rdMD.CalcExactMolWt(mol)
    self.assertTrue(feq(amw, 16.031, .001))

  def testPairValues(self):
    import base64
    testD = (
      ('CCCO',
       b'AQAAAAQAAAAAAIAABgAAACGECAABAAAAIoQIAAEAAABBhAgAAQAAACNEGAABAAAAQUQYAAEAAABC\nRBgAAQAAAA==\n'
       ),
      ('CNc1ccco1',
       b'AQAAAAQAAAAAAIAAEAAAACOECgABAAAAJIQKAAIAAABBhQoAAgAAAEKFCgABAAAAIsQKAAEAAABB\nxQoAAQAAAELFCgACAAAAIYQQAAEAAABChRAAAQAAAEOFEAACAAAAYYUQAAEAAAAjhBoAAQAAAEGF\nGgABAAAAQoUaAAIAAABhhRoAAQAAAEKIGgABAAAA\n'
       ),
    )
    for smi, txt in testD:
      pkl = base64.decodestring(txt)
      fp = rdMD.GetAtomPairFingerprint(Chem.MolFromSmiles(smi))
      fp2 = DataStructs.IntSparseIntVect(pkl)
      self.assertEqual(DataStructs.DiceSimilarity(fp, fp2), 1.0)
      self.assertEqual(fp, fp2)

  def testTorsionValues(self):
    import base64
    testD = (
      ('CCCO', b'AQAAAAgAAAD/////DwAAAAEAAAAAAAAAIECAAAMAAAABAAAA\n'),
      ('CNc1ccco1',
       b'AQAAAAgAAAD/////DwAAAAkAAAAAAAAAIICkSAEAAAABAAAAKVKgSQEAAAABAAAAKVCgUAEAAAAB\nAAAAKVCgUQEAAAABAAAAKVCkCAIAAAABAAAAKdCkCAIAAAABAAAAKVCgSAMAAAABAAAAKVCkSAMA\nAAABAAAAIICkSAMAAAABAAAA\n'
       ),
    )
    for smi, txt in testD:
      pkl = base64.decodestring(txt)
      fp = rdMD.GetTopologicalTorsionFingerprint(Chem.MolFromSmiles(smi))
      fp2 = DataStructs.LongSparseIntVect(pkl)
      self.assertEqual(DataStructs.DiceSimilarity(fp, fp2), 1.0)
      self.assertEqual(fp, fp2)

  def testAtomPairOptions(self):
    m1 = Chem.MolFromSmiles('c1ccccc1')
    m2 = Chem.MolFromSmiles('c1ccccn1')

    fp1 = rdMD.GetAtomPairFingerprint(m1)
    fp2 = rdMD.GetAtomPairFingerprint(m2)
    self.assertNotEqual(fp1, fp2)

    fp1 = rdMD.GetAtomPairFingerprint(m1, atomInvariants=[1] * 6)
    fp2 = rdMD.GetAtomPairFingerprint(m2, atomInvariants=[1] * 6)
    self.assertEqual(fp1, fp2)

    fp1 = rdMD.GetAtomPairFingerprint(m1, atomInvariants=[1] * 6)
    fp2 = rdMD.GetAtomPairFingerprint(m2, atomInvariants=[2] * 6)
    self.assertNotEqual(fp1, fp2)

    fp1 = rdMD.GetHashedAtomPairFingerprintAsBitVect(m1)
    fp2 = rdMD.GetHashedAtomPairFingerprintAsBitVect(m2)
    self.assertNotEqual(fp1, fp2)

    fp1 = rdMD.GetHashedAtomPairFingerprintAsBitVect(m1, atomInvariants=[1] * 6)
    fp2 = rdMD.GetHashedAtomPairFingerprintAsBitVect(m2, atomInvariants=[1] * 6)
    self.assertEqual(fp1, fp2)

    fp1 = rdMD.GetHashedAtomPairFingerprintAsBitVect(m1, atomInvariants=[1] * 6)
    fp2 = rdMD.GetHashedAtomPairFingerprintAsBitVect(m2, atomInvariants=[2] * 6)
    self.assertNotEqual(fp1, fp2)

    fp1 = rdMD.GetTopologicalTorsionFingerprint(m1)
    fp2 = rdMD.GetTopologicalTorsionFingerprint(m2)
    self.assertNotEqual(fp1, fp2)

    fp1 = rdMD.GetTopologicalTorsionFingerprint(m1, atomInvariants=[1] * 6)
    fp2 = rdMD.GetTopologicalTorsionFingerprint(m2, atomInvariants=[1] * 6)
    self.assertEqual(fp1, fp2)

    fp1 = rdMD.GetTopologicalTorsionFingerprint(m1, atomInvariants=[1] * 6)
    fp2 = rdMD.GetTopologicalTorsionFingerprint(m2, atomInvariants=[2] * 6)
    self.assertNotEqual(fp1, fp2)

    fp1 = rdMD.GetHashedTopologicalTorsionFingerprintAsBitVect(m1)
    fp2 = rdMD.GetHashedTopologicalTorsionFingerprintAsBitVect(m2)
    self.assertNotEqual(fp1, fp2)

    fp1 = rdMD.GetHashedTopologicalTorsionFingerprintAsBitVect(m1, atomInvariants=[1] * 6)
    fp2 = rdMD.GetHashedTopologicalTorsionFingerprintAsBitVect(m2, atomInvariants=[1] * 6)
    self.assertEqual(fp1, fp2)

    fp1 = rdMD.GetHashedTopologicalTorsionFingerprintAsBitVect(m1, atomInvariants=[1] * 6)
    fp2 = rdMD.GetHashedTopologicalTorsionFingerprintAsBitVect(m2, atomInvariants=[2] * 6)
    self.assertNotEqual(fp1, fp2)

  def testMolFormula(self):
    m = Chem.MolFromSmiles("[2H]C([3H])O")
    formula = rdMD.CalcMolFormula(m)
    self.assertEqual(formula, 'CH4O')
    formula = rdMD.CalcMolFormula(m, separateIsotopes=True)
    self.assertEqual(formula, 'CH2DTO')
    formula = rdMD.CalcMolFormula(m, separateIsotopes=True, abbreviateHIsotopes=False)
    self.assertEqual(formula, 'CH2[2H][3H]O')

    m = Chem.MolFromSmiles("[2H][13CH2]CO")
    formula = rdMD.CalcMolFormula(m)
    self.assertEqual(formula, 'C2H6O')
    formula = rdMD.CalcMolFormula(m, separateIsotopes=True)
    self.assertEqual(formula, 'C[13C]H5DO')

  def testSpiroAndBridgeheads(self):
    m = Chem.MolFromSmiles("C1CC2CCC1CC2")
    self.assertEqual(rdMD.CalcNumSpiroAtoms(m), 0)
    sa = []
    self.assertEqual(rdMD.CalcNumSpiroAtoms(m, atoms=sa), 0)
    self.assertEqual(len(sa), 0)

    self.assertEqual(rdMD.CalcNumBridgeheadAtoms(m), 2)
    sa = []
    self.assertEqual(rdMD.CalcNumBridgeheadAtoms(m, atoms=sa), 2)
    self.assertEqual(len(sa), 2)
    self.assertEqual(sorted(sa), [2, 5])

    m = Chem.MolFromSmiles("C1CCC2(C1)CC1CCC2CC1")
    self.assertEqual(rdMD.CalcNumSpiroAtoms(m), 1)
    sa = []
    self.assertEqual(rdMD.CalcNumSpiroAtoms(m, atoms=sa), 1)
    self.assertEqual(len(sa), 1)
    self.assertEqual(sorted(sa), [3])

    self.assertEqual(rdMD.CalcNumBridgeheadAtoms(m), 2)
    sa = []
    self.assertEqual(rdMD.CalcNumBridgeheadAtoms(m, atoms=sa), 2)
    self.assertEqual(len(sa), 2)
    self.assertEqual(sorted(sa), [6, 9])

  def testNumRotatableBonds(self):
    for s in [
        "C1CC1CC",
        "CCNC(=O)NCC",
        'Cc1cccc(C)c1c1c(C)cccc1C',
        'CCc1cccc(C)c1c1c(C)cccc1CC',
        'Cc1cccc(C)c1c1c(C)nccc1C',
        'Cc1cccc(C)c1c1c(C)cccc1',
        'CCO',
    ]:

      m = Chem.MolFromSmiles(s)

      v1 = rdMD.CalcNumRotatableBonds(m)

      v2 = rdMD.CalcNumRotatableBonds(m, False)
      v3 = rdMD.CalcNumRotatableBonds(m, True)

      v4 = rdMD.CalcNumRotatableBonds(m, rdMD.NumRotatableBondsOptions.Default)
      v5 = rdMD.CalcNumRotatableBonds(m, rdMD.NumRotatableBondsOptions.NonStrict)
      v6 = rdMD.CalcNumRotatableBonds(m, rdMD.NumRotatableBondsOptions.Strict)
      v7 = rdMD.CalcNumRotatableBonds(m, rdMD.NumRotatableBondsOptions.StrictLinkages)

      self.assertEqual(v1, v4)
      self.assertEqual(v2, v5)
      self.assertEqual(v3, v6)

  def testProperties(self):
    props = rdMD.Properties()
    names = list(props.GetAvailableProperties())
    self.assertEqual(names, list(props.GetPropertyNames()))
    m = Chem.MolFromSmiles("C1CC1CC")
    results = props.ComputeProperties(m)

    for i, name in enumerate(names):
      props = rdMD.Properties([name])
      res = props.ComputeProperties(m)
      self.assertEqual(len(res), 1)
      self.assertEqual(res[0], results[i])
      self.assertEqual(props.GetPropertyNames()[0], names[i])
      self.assertEqual(len(props.GetPropertyNames()), 1)

    try:
      props = rdMD.Properties([1, 2, 3])
      self.assertEqual("should not get here", "but did")
    except TypeError:
      pass

    try:
      props = rdMD.Properties(["property that doesn't exist"])
      self.assertEqual("should not get here", "but did")
    except RuntimeError:
      pass

  def testPythonDescriptorFunctor(self):

    class NumAtoms(Descriptors.PropertyFunctor):

      def __init__(self):
        Descriptors.PropertyFunctor.__init__(self, "NumAtoms", "1.0.0")

      def __call__(self, mol):
        return mol.GetNumAtoms()

    numAtoms = NumAtoms()
    rdMD.Properties.RegisterProperty(numAtoms)
    props = rdMD.Properties(["NumAtoms"])
    self.assertEqual(1, props.ComputeProperties(Chem.MolFromSmiles("C"))[0])

    self.assertTrue("NumAtoms" in rdMD.Properties.GetAvailableProperties())
    # check memory
    del numAtoms
    self.assertEqual(1, props.ComputeProperties(Chem.MolFromSmiles("C"))[0])
    self.assertTrue("NumAtoms" in rdMD.Properties.GetAvailableProperties())

    m = Chem.MolFromSmiles("c1ccccc1")
    properties = rdMD.Properties()
    for name, value in zip(properties.GetPropertyNames(), properties.ComputeProperties(m)):
      print(name, value)

    properties = rdMD.Properties(['exactmw', 'lipinskiHBA'])
    for name, value in zip(properties.GetPropertyNames(), properties.ComputeProperties(m)):
      print(name, value)

  def testPropertyRanges(self):
    query = rdMD.MakePropertyRangeQuery("exactmw", 0, 1000)
    self.assertTrue(query.Match(Chem.MolFromSmiles("C")))

    query = rdMD.MakePropertyRangeQuery("exactmw", 1000, 10000)
    self.assertFalse(query.Match(Chem.MolFromSmiles("C")))

  def testNumStereoCenters(self):
    m = Chem.MolFromSmiles('CC(F)(Cl)[C@H](Cl)Br')
    self.assertEqual(rdMD.CalcNumAtomStereoCenters(m), 2)
    self.assertEqual(rdMD.CalcNumUnspecifiedAtomStereoCenters(m), 1)
    # Tests from Berend Huisman:
    for (smiles, expected) in (
      ("C", 0),
      ("c1ccccc1", 0),
      ("CC(Cl)Br", 1),
      ("CCC(C)C(Cl)Br", 2),
      ("CCC(C(Cl)Br)C(F)I", 3),
      ("[H][C@](F)(I)C(CC)C(Cl)Br", 3),
      ("[H][C@](F)(I)[C@@]([H])(CC)C(Cl)Br", 3),
    ):
      mol = Chem.MolFromSmiles(smiles)
      actual = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
      self.assertEqual(rdMD.CalcNumAtomStereoCenters(mol), expected)
    for (smiles, expected) in (
      ("C", 0),
      ("c1ccccc1", 0),
      ("CC(Cl)Br", 1),
      ("CCC(C)C(Cl)Br", 2),
      ("CCC(C(Cl)Br)C(F)I", 3),
      ("[H][C@](F)(I)C(CC)C(Cl)Br", 2),
      ("[H][C@](F)(I)[C@@]([H])(CC)C(Cl)Br", 1),
    ):
      mol = Chem.MolFromSmiles(smiles)
      actual = sum(1 for x in Chem.FindMolChiralCenters(mol, includeUnassigned=True) if x[1] == '?')
      self.assertEqual(actual, expected)
      self.assertEqual(rdMD.CalcNumUnspecifiedAtomStereoCenters(mol), expected)

  def testGithub1749(self):
    mol = Chem.MolFromSmiles("c1ccccc1O")
    self.assertRaises(ValueError,
                      lambda: rdMD.GetMorganFingerprintAsBitVect(mol, 2, fromAtoms=[10]))

  def testCustomVSA(self):
    mol = Chem.MolFromSmiles("c1ccccc1O")
    peoe_vsa = rdMD.PEOE_VSA_(mol)
    AllChem.ComputeGasteigerCharges(mol)
    bins = [-.3, -.25, -.20, -.15, -.10, -.05, 0, .05, .10, .15, .20, .25, .30]
    custom_vsa = rdMD.CustomProp_VSA_(mol, customPropName='_GasteigerCharge', bins=bins)
    for p, c in zip(peoe_vsa, custom_vsa):
      self.assertTrue(feq(p, c, .001))

  def testGithub1973(self):

    smiles = ("c1ccccc1S", "c1cscc1", "CC(=S)C", "CSC", "CS(=O)C", "CP(C)C", "CP=O", "CP(C)(C)=O",
              "C[PH](C)=O")
    orig_tpsa = (0, 0, 0, 0, 17.07, 0.0, 17.07, 17.07, 17.07)
    new_tpsa = (38.8, 28.24, 32.09, 25.30, 36.28, 13.59, 51.21, 26.88, 40.54)
    for i, smi in enumerate(smiles):
      mol = Chem.MolFromSmiles(smi)
      oTPSA = rdMD.CalcTPSA(mol)
      self.assertAlmostEqual(oTPSA, orig_tpsa[i], 2)
      nTPSA = rdMD.CalcTPSA(mol, force=True, includeSandP=True)
      self.assertAlmostEqual(nTPSA, new_tpsa[i], 2)


if __name__ == '__main__':
  unittest.main()
