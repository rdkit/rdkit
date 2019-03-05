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
from rdkit.Chem import rdDepictor
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
    res = np.zeros(nl, np.float)

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


if __name__ == '__main__':
    rdDepictor.SetPreferCoordGen(False)
    unittest.main()
