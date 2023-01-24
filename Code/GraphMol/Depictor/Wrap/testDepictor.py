# Automatically adapted for numpy.oldnumeric Jun 27, 2008 by -c

#
#  $Id: testDepictor.py 2112 2012-07-02 09:47:45Z glandrum $
#
# pylint:disable=E1101,C0111,C0103,R0904
import unittest
import os
import sys
import numpy as np

from rdkit import Chem
from rdkit.Chem import rdDepictor, rdMolAlign
from rdkit import Geometry
from rdkit import RDConfig
from rdkit.Chem.ChemUtils import AlignDepict


def feq(v1, v2, tol2=1e-4):
    return abs(v1 - v2) <= tol2


def ptEq(pt1, pt2, tol=1e-4):
    return feq(pt1.x, pt2.x, tol) and feq(pt1.y, pt2.y, tol) and feq(pt1.z, pt2.z, tol)


def getDistMat(mol):
    conf = mol.GetConformer()
    nat = mol.GetNumAtoms()
    nl = nat * (nat - 1) // 2
    res = np.zeros(nl, float)

    for i in range(1, nat):
        pi = conf.GetAtomPosition(i)
        idx = i * (i - 1) // 2
        for j in range(i):
            pj = conf.GetAtomPosition(j)
            pj -= pi
            res[idx + j] = pj.Length()

    return res


def compareCoords(m, molFile):
    mo = Chem.MolFromMolFile(molFile)
    co = mo.GetConformer()

    ci = m.GetConformer()
    nat = m.GetNumAtoms()
    if (nat != mo.GetNumAtoms()):
        return 0

    for i in range(nat):
        pos = ci.GetAtomPosition(i)
        opos = co.GetAtomPosition(i)
        if not ptEq(pos, opos):
            print(Chem.MolToMolBlock(m))
            print(Chem.MolToMolBlock(mo))
            return 0
    return 1


def compareWithOld(smilesFile, sdFile):
    smiSup = Chem.SmilesMolSupplier(smilesFile, ",", 0, -1)
    sdsup = Chem.SDMolSupplier(sdFile)
    im = 0
    for mol in smiSup:
        omol = sdsup[im]
        rdDepictor.Compute2DCoords(mol, canonOrient=False)
        conf = mol.GetConformer()
        oconf = omol.GetConformer()
        nat = mol.GetNumAtoms()
        for i in range(nat):
            pos = conf.GetAtomPosition(i)
            opos = oconf.GetAtomPosition(i)
            if not ptEq(pos, opos):
                print(Chem.MolToMolBlock(omol), file=sys.stderr)
                print('> <Failed>\n%d\n' % i, file=sys.stderr)
                print("$$$$", file=sys.stderr)
                print(Chem.MolToMolBlock(mol), file=sys.stderr)
                print('> <Failed>\n%d\n' % i, file=sys.stderr)
                print("$$$$", file=sys.stderr)
                return 0
        im += 1
    return 1


def stereoCompare(smilesFile):
    smiSup = Chem.SmilesMolSupplier(smilesFile, ",", 0, -1)
    for mol in smiSup:
        rdDepictor.Compute2DCoords(mol, canonOrient=False)
        mb = Chem.MolToMolBlock(mol)
        nmol = Chem.MolFromMolBlock(mb)
        matches = nmol.GetSubstructMatches(mol, False)
        dbnds = [x for x in mol.GetBonds() if (x.GetBondType() == Chem.BondType.DOUBLE and
                                               x.GetStereo() > Chem.BondStereo.STEREOANY) ]
        ok = True
        for match in matches:
            for bnd in dbnds:
                obnd = nmol.GetBondBetweenAtoms(
                    match[bnd.GetBeginAtomIdx()], match[bnd.GetEndAtomIdx()])
                assert (obnd.GetBondType() == Chem.BondType.DOUBLE)
            if ok:
                break
        if not ok:
            print(Chem.MolToMolBlock(mol), file=sys.stderr)
            print("$$$$", file=sys.stderr)
            return 0
    return 1


class TestCase(unittest.TestCase):

    def _test0First200(self):
        # this test is disabled because it's not particularly useful and
        # causes problems every time anything changes.
        fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'Depictor', 'test_data',
                             'first_200.tpsa.csv')
        #smiSup = Chem.SmilesMolSupplier(fileN, ",", 0, -1)

        ofile = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'Depictor', 'test_data',
                             'first_200.python.sdf')
        self.assertTrue(compareWithOld(fileN, ofile))

    def test1CisTrans(self):
        fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'Depictor', 'test_data',
                             "cis_trans_cases.csv")
        self.assertTrue(stereoCompare(fileN))

    def test2Coords(self):
        m1 = Chem.MolFromSmiles('C1CCC1CC')
        coordMap = {0: Geometry.Point2D(0, 0),
                    1: Geometry.Point2D(1.5, 0),
                    2: Geometry.Point2D(1.5, 1.5),
                    3: Geometry.Point2D(0, 1.5)}
        rdDepictor.Compute2DCoords(m1, coordMap=coordMap)
        conf = m1.GetConformer(0)
        for i in range(4):
            self.assertTrue(
              ptEq(conf.GetAtomPosition(i), Geometry.Point3D(coordMap[i].x, coordMap[i].y, 0.0)))

        m1 = Chem.MolFromSmiles('CCC')
        try:
            rdDepictor.Compute2DCoords(m1, coordMap=coordMap)
            ok = 0
        except ValueError:
            ok = 1
        self.assertTrue(ok)

    def test3IssueSF1526844(self):
        t = Chem.MolFromSmiles('c1nc(N)ccc1')
        rdDepictor.Compute2DCoords(t, canonOrient=False)

        m2 = Chem.MolFromSmiles('c1nc(NC=O)ccc1')
        AlignDepict.AlignDepict(m2, t)
        expected = [Geometry.Point3D(1.5, 0.0, 0.0), Geometry.Point3D(0.75, -1.299, 0.0),
                    Geometry.Point3D(-0.75, -1.299, 0.0), Geometry.Point3D(-1.5, -2.5981, 0.0),
                    Geometry.Point3D(-3.0, -2.5981, 0.0), Geometry.Point3D(-3.75, -3.8971, 0.0),
                    Geometry.Point3D(-1.5, 0.0, 0.0), Geometry.Point3D(-0.75, 1.2990, 0.0),
                    Geometry.Point3D(0.75, 1.2990, 0.0)]

        nat = m2.GetNumAtoms()
        conf = m2.GetConformer()
        for i in range(nat):
            pos = conf.GetAtomPosition(i)
            self.assertTrue(ptEq(pos, expected[i], 0.001))

    def test4SamplingSpread(self):
        mol = Chem.MolFromMolFile(
          os.path.join(RDConfig.RDBaseDir, 'Code/GraphMol/Depictor', 'test_data/7UPJ_xtal.mol'))

        # default mode
        rdDepictor.Compute2DCoords(mol, canonOrient=False)
        self.assertTrue(
          compareCoords(mol, os.path.join(RDConfig.RDBaseDir, 'Code/GraphMol/Depictor',
                                          'test_data/7UPJ_default.mol')))

        # spread the structure as much as possible by sampling
        rdDepictor.Compute2DCoords(mol, canonOrient=False, nFlipsPerSample=3, nSample=100,
                                   sampleSeed=100, permuteDeg4Nodes=1)
        self.assertTrue(
          compareCoords(mol, os.path.join(RDConfig.RDBaseDir, 'Code/GraphMol/Depictor',
                                          'test_data/7UPJ_spread.mol')))

    def test5SamplingMimic3D(self):
        mol = Chem.MolFromMolFile(
          os.path.join(RDConfig.RDBaseDir, 'Code/GraphMol/Depictor', 'test_data/7UPJ_xtal.mol'))
        dmat3D = getDistMat(mol)

        # now mimic the coordinate with a very small weight
        rdDepictor.Compute2DCoordsMimicDistmat(mol, dmat3D, weightDistMat=0.001)
        self.assertTrue(
          compareCoords(mol, os.path.join(RDConfig.RDBaseDir, 'Code/GraphMol/Depictor',
                                          'test_data/7UPJ_mimic3D_1.mol')))

        # now mimic the coordinate with a very small weight
        rdDepictor.Compute2DCoordsMimicDistmat(mol, dmat3D, weightDistMat=0.003)
        self.assertTrue(
          compareCoords(mol, os.path.join(RDConfig.RDBaseDir, 'Code/GraphMol/Depictor',
                                          'test_data/7UPJ_mimic3D_2.mol')))

        #mb = Chem.MolToMolBlock(mol)
        #ofile = open('../test_data/7UPJ_mimic3D_2.mol', 'w')
        # ofile.write(mb)
        # ofile.close()

    def test6ChangeBondLength(self):
        m = Chem.MolFromSmiles('CC')
        rdDepictor.Compute2DCoords(m)
        conf = m.GetConformer()
        self.assertAlmostEqual(conf.GetAtomPosition(0).x, -0.750, 3)
        self.assertAlmostEqual(conf.GetAtomPosition(1).x, 0.750, 3)
        rdDepictor.Compute2DCoords(m, bondLength=1.0)
        conf = m.GetConformer()
        self.assertAlmostEqual(conf.GetAtomPosition(0).x, -0.500, 3)
        self.assertAlmostEqual(conf.GetAtomPosition(1).x, 0.500, 3)
        rdDepictor.Compute2DCoords(m)
        conf = m.GetConformer()
        self.assertAlmostEqual(conf.GetAtomPosition(0).x, -0.750, 3)
        self.assertAlmostEqual(conf.GetAtomPosition(1).x, 0.750, 3)

    def testConstrainedCoords(self):
        templ = Chem.MolFromSmiles('c1nccc2n1ccc2')
        rdDepictor.Compute2DCoords(templ)
        m1 = Chem.MolFromSmiles('c1cccc2ncn3cccc3c21')
        rdDepictor.GenerateDepictionMatching2DStructure(m1, templ)
        m2 = Chem.MolFromSmiles('c1cc(Cl)cc2ncn3cccc3c21')
        rdDepictor.Compute2DCoords(m2)
        refPatt1 = Chem.MolFromSmarts('*1****2*1***2')
        rdDepictor.GenerateDepictionMatching2DStructure(m2, templ, -1, refPatt1)
        fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'Depictor', 'test_data',
                             '1XP0_ligand.sdf')

        xp0_lig = Chem.MolFromMolFile(fileN)
        xp0_lig_2d = Chem.Mol(xp0_lig)
        rdDepictor.GenerateDepictionMatching3DStructure(xp0_lig_2d, xp0_lig)
        xp0_ref = Chem.MolFromSmarts('[#6]1~[#7][#6]~[#6]2[#6](=[#8])[#7]~[#6](c3ccccc3)[#7][#7]12')
        rdDepictor.GenerateDepictionMatching3DStructure(xp0_lig_2d, xp0_lig, -1, xp0_ref)


    def testGenerate2DDepictionRefPatternAtomMap(self):
        indazoleMolblock = """
     RDKit          2D

  9 10  0  0  0  0  0  0  0  0999 V2000
   -6.0878    2.4335    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3867    1.6835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -7.3867    0.1833    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -6.0878   -0.5666    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7887    0.1833    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -4.7887    1.6835    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4897   -0.5664    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1906    1.6833    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1906    0.1835    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  8  9  2  0
  6  8  1  0
  7  9  1  0
  7  5  1  0
M  END"""
        indazoleRef = Chem.MolFromMolBlock(indazoleMolblock)
        cycloheptylPyrazole = Chem.MolFromSmiles("c1cc(C2CCCCCC2)[nH]n1")

        # test using refPattern
        refPatt = Chem.MolFromSmarts("a1aan[nH]1")
        rdDepictor.GenerateDepictionMatching2DStructure(cycloheptylPyrazole, indazoleRef, refPatt=refPatt)
        self.assertEqual(cycloheptylPyrazole.GetNumConformers(), 1)
        molMatchVect = cycloheptylPyrazole.GetSubstructMatch(refPatt)
        self.assertEqual(len(molMatchVect), refPatt.GetNumAtoms())
        refMatchVect = indazoleRef.GetSubstructMatch(refPatt)
        self.assertEqual(len(refMatchVect), refPatt.GetNumAtoms())
        atomMap = tuple(zip(refMatchVect, molMatchVect))
        msd = 0.0
        for refIdx, molIdx in atomMap:
            msd += (indazoleRef.GetConformer().GetAtomPosition(refIdx) -
                    cycloheptylPyrazole.GetConformer().GetAtomPosition(molIdx)).LengthSq()
        msd /= len(molMatchVect)
        self.assertAlmostEqual(msd, 0.0)
        # try with a pattern larger than the reference molecule
        hugePatt = Chem.MolFromSmarts("CCCCCCCCCCCCCCCCCCCCCCCCCCC")
        with self.assertRaises(ValueError):
            rdDepictor.GenerateDepictionMatching2DStructure(
                cycloheptylPyrazole, indazoleRef, refPatt=hugePatt)

        # try with an out of range confId
        with self.assertRaises(ValueError):
            rdDepictor.GenerateDepictionMatching2DStructure(
                cycloheptylPyrazole, indazoleRef, confId=1, refPatt=refPatt)

        # test using atomMap directly
        cycloheptylPyrazole.RemoveAllConformers()
        rdDepictor.GenerateDepictionMatching2DStructure(cycloheptylPyrazole, indazoleRef, atomMap=atomMap)
        self.assertEqual(cycloheptylPyrazole.GetNumConformers(), 1)
        msd = 0.0
        for refIdx, molIdx in atomMap:
            msd += (indazoleRef.GetConformer().GetAtomPosition(refIdx) -
                    cycloheptylPyrazole.GetConformer().GetAtomPosition(molIdx)).LengthSq()
        msd /= len(atomMap)
        self.assertAlmostEqual(msd, 0.0)

        # try with an atomMap larger than the reference molecule
        atomMapHuge = list(atomMap) + [(0, 0) for i in range(indazoleRef.GetNumAtoms())]
        with self.assertRaises(ValueError):
            rdDepictor.GenerateDepictionMatching2DStructure(
                cycloheptylPyrazole, indazoleRef, atomMap=atomMapHuge)

        # try with an atomMap with out of range indices
        atomMapOutOfRange = list(atomMap) + [(100, 100)]
        with self.assertRaises(ValueError):
            rdDepictor.GenerateDepictionMatching2DStructure(
                cycloheptylPyrazole, indazoleRef, atomMap=atomMapOutOfRange)

        # try with an out of range confId
        with self.assertRaises(ValueError):
            rdDepictor.GenerateDepictionMatching2DStructure(
                cycloheptylPyrazole, indazoleRef, atomMap=atomMap, confId=1)

    def testGenerate2DDepictionAllowRGroups(self):
        templateMolblock = """
     RDKit          2D

  9  9  0  0  0  0  0  0  0  0999 V2000
   -0.8929    1.0942    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1919    0.3442    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1919   -1.1558    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8929   -1.9059    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4060   -1.1558    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4060    0.3442    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4910    1.0942    0.0000 R1  0  0  0  0  0  0  0  0  0  0  0  0
    1.7051    1.0942    0.0000 R2  0  0  0  0  0  0  0  0  0  0  0  0
   -3.4910   -1.9059    0.0000 R3  0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
  6  8  1  0
  3  9  1  0
  2  7  1  0
M  RGP  3   7   1   8   2   9   3
M  END"""
        templateRef = Chem.MolFromMolBlock(templateMolblock)
        orthoMeta = Chem.MolFromSmiles("c1ccc(-c2ccc(-c3ccccc3)c(-c3ccccc3)c2)cc1")
        ortho = Chem.MolFromSmiles("c1ccc(-c2ccccc2-c2ccccc2)cc1")
        meta = Chem.MolFromSmiles("c1ccc(-c2cccc(-c3ccccc3)c2)cc1")
        biphenyl = Chem.MolFromSmiles("c1ccccc1-c1ccccc1")
        phenyl = Chem.MolFromSmiles("c1ccccc1")

        atomMap = rdDepictor.GenerateDepictionMatching2DStructure(orthoMeta, templateRef)
        self.assertEqual(orthoMeta.GetNumConformers(), 1)

        for mol in (ortho, meta, biphenyl, phenyl):
            # fails as does not match template
            with self.assertRaises(ValueError):
                rdDepictor.GenerateDepictionMatching2DStructure(mol, templateRef)

            # succeeds with allowRGroups=true
            atomMap = rdDepictor.GenerateDepictionMatching2DStructure(mol, templateRef, allowRGroups=True)
            self.assertEqual(mol.GetNumConformers(), 1)
            msd = 0.0
            for refIdx, molIdx in atomMap:
                msd += (templateRef.GetConformer().GetAtomPosition(refIdx) -
                        mol.GetConformer().GetAtomPosition(molIdx)).LengthSq()
            msd /= len(atomMap)
            self.assertAlmostEqual(msd, 0.0)

        # test that using a refPattern with R groups and a reference without works
        pyridineRef = Chem.MolFromMolBlock("""
     RDKit          2D

  6  6  0  0  0  0  0  0  0  0999 V2000
   -0.8929    1.0942    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1919    0.3442    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.1919   -1.1558    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8929   -1.9059    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4060   -1.1558    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.4060    0.3442    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  2  0
  2  3  1  0
  3  4  2  0
  4  5  1  0
  5  6  2  0
  6  1  1  0
M  END""")
        genericRefPatternWithRGroups = Chem.MolFromSmarts("[*:3]a1a([*:1])aa([*:2])aa1")

        for mol in (ortho, meta, biphenyl, phenyl):
            atomMap = rdDepictor.GenerateDepictionMatching2DStructure(
                mol, pyridineRef, refPatt=genericRefPatternWithRGroups, allowRGroups=True)
            self.assertEqual(mol.GetNumConformers(), 1)
            msd = 0.0
            for refIdx, molIdx in atomMap:
                msd += (pyridineRef.GetConformer().GetAtomPosition(refIdx) -
                        mol.GetConformer().GetAtomPosition(molIdx)).LengthSq()
            msd /= len(atomMap)
            self.assertAlmostEqual(msd, 0.0)

    def testNormalizeStraighten(self):
        noradrenalineMJ = Chem.MolFromMolBlock("""
  MJ201100                      

 12 12  0  0  1  0  0  0  0  0999 V2000
    2.2687    1.0716    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.4437    1.0716    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0312    0.3572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4437   -0.3572    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.2062    0.3572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2062   -0.3572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0312   -0.3572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4437   -1.0716    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.4437    0.3572    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.2687    0.3572    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0312    1.0716    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2062    1.0716    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  3  2  1  0  0  0  0
  3  4  1  6  0  0  0
  3  5  1  0  0  0  0
  5  6  2  0  0  0  0
  6  7  1  0  0  0  0
  7  8  1  0  0  0  0
  7  9  2  0  0  0  0
  9 10  1  0  0  0  0
  9 11  1  0  0  0  0
 11 12  2  0  0  0  0
  5 12  1  0  0  0  0
M  END)""")

        noradrenalineMJCopy = Chem.Mol(noradrenalineMJ)
        conformer0 = noradrenalineMJCopy.GetConformer(0)
        conformer1 = Chem.Conformer(conformer0)
        noradrenalineMJCopy.AddConformer(conformer1, True)
        conformer1 = noradrenalineMJCopy.GetConformer(1)
        self.assertLess(rdMolAlign.CalcRMS(noradrenalineMJ, noradrenalineMJCopy, 0, 0), 1.e-5)
        self.assertLess(rdMolAlign.CalcRMS(noradrenalineMJ, noradrenalineMJCopy, 0, 1), 1.e-5)
        scalingFactor = rdDepictor.NormalizeDepiction(noradrenalineMJCopy, 1)
        self.assertLess(rdMolAlign.CalcRMS(noradrenalineMJ, noradrenalineMJCopy, 0, 0), 1.e-5)
        self.assertGreater(rdMolAlign.CalcRMS(noradrenalineMJ, noradrenalineMJCopy, 0, 1), 1.e-5)
        self.assertAlmostEqual(scalingFactor, 1.875, 3)
        conformer2 = Chem.Conformer(conformer1)
        noradrenalineMJCopy.AddConformer(conformer2, True)
        conformer2 = noradrenalineMJCopy.GetConformer(2)
        bond10_11Conf0 = conformer0.GetAtomPosition(11) - conformer0.GetAtomPosition(10)
        self.assertAlmostEqual(bond10_11Conf0.x, 0.825, 3)
        self.assertAlmostEqual(bond10_11Conf0.y, 0., 3)
        bond10_11Conf1 = conformer1.GetAtomPosition(11) - conformer1.GetAtomPosition(10)
        self.assertAlmostEqual(bond10_11Conf1.x, 1.513, 3)
        self.assertAlmostEqual(bond10_11Conf1.y, -0.321, 3)
        rdDepictor.StraightenDepiction(noradrenalineMJCopy, 1)
        bond10_11Conf1 = conformer1.GetAtomPosition(11) - conformer1.GetAtomPosition(10)
        self.assertAlmostEqual(bond10_11Conf1.x, 1.340, 3)
        self.assertAlmostEqual(bond10_11Conf1.y, -0.773, 3)
        bond4_11Conf1 = conformer1.GetAtomPosition(11) - conformer1.GetAtomPosition(4)
        self.assertAlmostEqual(bond4_11Conf1.x, 0., 3)
        self.assertAlmostEqual(bond4_11Conf1.y, 1.547, 3)
        rdDepictor.StraightenDepiction(noradrenalineMJCopy, 2, True)
        bond10_11Conf2 = conformer2.GetAtomPosition(11) - conformer2.GetAtomPosition(10)
        self.assertAlmostEqual(bond10_11Conf2.x, 1.547, 3)
        self.assertAlmostEqual(bond10_11Conf2.y, 0.0, 3)
        bond4_11Conf2 = conformer2.GetAtomPosition(11) - conformer2.GetAtomPosition(4)
        self.assertAlmostEqual(bond4_11Conf2.x, -0.773, 3)
        self.assertAlmostEqual(bond4_11Conf2.y, 1.339, 3)

        noradrenalineMJCopy = Chem.Mol(noradrenalineMJ)
        conformer0 = noradrenalineMJCopy.GetConformer(0)
        conformer1 = Chem.Conformer(conformer0)
        noradrenalineMJCopy.AddConformer(conformer1, True)
        conformer1 = noradrenalineMJCopy.GetConformer(1)
        scalingFactor = rdDepictor.NormalizeDepiction(noradrenalineMJCopy, 1, -1)
        self.assertLess(rdMolAlign.CalcRMS(noradrenalineMJ, noradrenalineMJCopy, 0, 0), 1.e-5)
        self.assertGreater(rdMolAlign.CalcRMS(noradrenalineMJ, noradrenalineMJCopy, 0, 1), 1.e-5)
        self.assertAlmostEqual(scalingFactor, 1.875, 3)
        conformer2 = Chem.Conformer(conformer1)
        noradrenalineMJCopy.AddConformer(conformer2, True)
        conformer2 = noradrenalineMJCopy.GetConformer(2)
        bond10_11Conf0 = conformer0.GetAtomPosition(11) - conformer0.GetAtomPosition(10)
        self.assertAlmostEqual(bond10_11Conf0.x, 0.825, 3)
        self.assertAlmostEqual(bond10_11Conf0.y, 0., 3)
        bond10_11Conf1 = conformer1.GetAtomPosition(11) - conformer1.GetAtomPosition(10)
        self.assertAlmostEqual(bond10_11Conf1.x, 0.321, 3)
        self.assertAlmostEqual(bond10_11Conf1.y, 1.513, 3)
        rdDepictor.StraightenDepiction(noradrenalineMJCopy, 1)
        bond10_11Conf1 = conformer1.GetAtomPosition(11) - conformer1.GetAtomPosition(10)
        self.assertAlmostEqual(bond10_11Conf1.x, 0.0, 3)
        self.assertAlmostEqual(bond10_11Conf1.y, 1.547, 3)
        rdDepictor.StraightenDepiction(noradrenalineMJCopy, 2, True)
        bond10_11Conf2 = conformer2.GetAtomPosition(11) - conformer2.GetAtomPosition(10)
        self.assertAlmostEqual(bond10_11Conf2.x, bond10_11Conf1.x, 3)
        self.assertAlmostEqual(bond10_11Conf2.y, bond10_11Conf1.y, 3)

        noradrenalineMJCopy = Chem.Mol(noradrenalineMJ)
        conformer0 = noradrenalineMJCopy.GetConformer(0)
        conformer1 = Chem.Conformer(conformer0)
        noradrenalineMJCopy.AddConformer(conformer1, True)
        conformer1 = noradrenalineMJCopy.GetConformer(1)
        scalingFactor = rdDepictor.NormalizeDepiction(noradrenalineMJCopy, 1, 0, 3.0)
        self.assertLess(rdMolAlign.CalcRMS(noradrenalineMJ, noradrenalineMJCopy, 0, 0), 1.e-5)
        self.assertGreater(rdMolAlign.CalcRMS(noradrenalineMJ, noradrenalineMJCopy, 0, 1), 1.e-5)
        self.assertAlmostEqual(scalingFactor, 3.0, 3)
        conformer2 = Chem.Conformer(conformer1)
        noradrenalineMJCopy.AddConformer(conformer2, True)
        conformer2 = noradrenalineMJCopy.GetConformer(2)
        conformer3 = Chem.Conformer(conformer1)
        noradrenalineMJCopy.AddConformer(conformer3, True)
        conformer3 = noradrenalineMJCopy.GetConformer(3)
        bond10_11Conf0 = conformer0.GetAtomPosition(11) - conformer0.GetAtomPosition(10)
        self.assertAlmostEqual(bond10_11Conf0.x, 0.825, 3)
        self.assertAlmostEqual(bond10_11Conf0.y, 0., 3)
        bond10_11Conf1 = conformer1.GetAtomPosition(11) - conformer1.GetAtomPosition(10)
        self.assertAlmostEqual(bond10_11Conf1.x, 2.475, 3)
        self.assertAlmostEqual(bond10_11Conf1.y, 0., 3)
        rdDepictor.StraightenDepiction(noradrenalineMJCopy, 1)
        bond10_11Conf1 = conformer1.GetAtomPosition(11) - conformer1.GetAtomPosition(10)
        self.assertAlmostEqual(bond10_11Conf1.x, 2.143, 3)
        self.assertAlmostEqual(bond10_11Conf1.y, -1.237, 3)
        bond4_11Conf1 = conformer1.GetAtomPosition(11) - conformer1.GetAtomPosition(4)
        self.assertAlmostEqual(bond4_11Conf1.x, 0., 3)
        self.assertAlmostEqual(bond4_11Conf1.y, 2.475, 3)
        rdDepictor.StraightenDepiction(noradrenalineMJCopy, 2, True)
        bond10_11Conf2 = conformer2.GetAtomPosition(11) - conformer2.GetAtomPosition(10)
        bond10_11Conf3 = conformer3.GetAtomPosition(11) - conformer3.GetAtomPosition(10)
        self.assertAlmostEqual(bond10_11Conf2.x, bond10_11Conf3.x, 3)
        self.assertAlmostEqual(bond10_11Conf2.y, bond10_11Conf3.y, 3)
        bond4_11Conf2 = conformer2.GetAtomPosition(11) - conformer2.GetAtomPosition(4)
        bond4_11Conf3 = conformer3.GetAtomPosition(11) - conformer3.GetAtomPosition(4)
        self.assertAlmostEqual(bond4_11Conf2.x, bond4_11Conf3.x, 3)
        self.assertAlmostEqual(bond4_11Conf2.y, bond4_11Conf3.y, 3)

if __name__ == '__main__':
    rdDepictor.SetPreferCoordGen(False)
    unittest.main()
