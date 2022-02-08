
import unittest
import os
import copy
import math
import numpy

import pickle

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDistGeom, ChemicalForceFields, rdMolAlign
import rdkit.DistanceGeometry as DG
from rdkit import RDConfig, rdBase
from rdkit.Geometry import rdGeometry as geom
from rdkit.Geometry import ComputeSignedDihedralAngle
from rdkit.RDLogger import logger
logger = logger()


def feq(v1, v2, tol=1.e-4):
    return abs(v1 - v2) < tol


def lstEq(l1, l2, tol=1.0e-4):
    ln = len(l1)
    if (ln != len(l2)):
        return 0

    for i in range(ln):
        if abs(l1[i] - l2[i]) > tol:
            return 0
    return 1


def compareWithOld(smilesFile, sdFile):
    smiSup = Chem.SmilesMolSupplier(smilesFile, ",", 0, -1)
    sdsup = Chem.SDMolSupplier(sdFile)
    im = 0
    for mol in smiSup:
        cid = rdDistGeom.EmbedMolecule(mol, 10, 1)
        omol = sdsup[im]
        assert cid == 0
        conf = mol.GetConformer(0)
        oconf = omol.GetConformer()
        nat = mol.GetNumAtoms()
        for i in range(nat):
            #atm = mol.GetAtomWithIdx(i)
            #oatm = omol.GetAtomWithIdx(i)
            pos = conf.GetAtomPosition(i)
            opos = oconf.GetAtomPosition(i)
            if not lstEq(pos, opos):
                return 0
        im += 1
    return 1


def compareMatrices(bm1, bm2, map, tol=1.0e-5):
    N = numpy.shape(bm1)[0]
    for i in range(1, N):
        for j in range(i):
            l, m = map[i], map[j]
            if (l < m):
                l, m = m, l
            if (abs(bm1[l, m] - bm2[i, j]) > tol):
                return 0

            if (abs(bm1[m, l] - bm2[j, i]) > tol):
                return 0

    return 1


def compareOrder(smi1, smi2, tol=1.0e-5):
    m1 = Chem.MolFromSmiles(smi1)
    m2 = Chem.MolFromSmiles(smi2)
    bm1 = rdDistGeom.GetMoleculeBoundsMatrix(m1)
    bm2 = rdDistGeom.GetMoleculeBoundsMatrix(m2)
    map = m1.GetSubstructMatch(m2)
    return compareMatrices(bm1, bm2, map, tol)


def computeDist(lst1, lst2):
    res = 0.0
    for i, val in enumerate(lst1):
        res += (val - lst2[i]) * (val - lst2[i])
    res = math.sqrt(res)
    return res


def computeChiralVol(pt1, pt2, pt3, pt4):
    v1 = pt1 - pt4
    v2 = pt2 - pt4
    v3 = pt3 - pt4
    cp = v2.CrossProduct(v3)
    vol = v1.DotProduct(cp)
    return vol


class TestCase(unittest.TestCase):

    def setUp(self):
        pass

    def _test0Cdk2(self):
        fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'DistGeomHelpers', 'test_data',
                             'cis_trans_cases.csv')

        ofile = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'DistGeomHelpers', 'test_data',
                             'embedDistOpti.sdf')
        self.assertTrue(compareWithOld(fileN, ofile))

    def test1Small(self):
        #writer = Chem.SDWriter("test.sdf")
        # single double and tripple atoms cases should not fail
        mol = Chem.MolFromSmiles('O')
        rdDistGeom.EmbedMolecule(mol, 10, 1)
        conf = mol.GetConformer()
        self.assertTrue(lstEq(conf.GetAtomPosition(0), [0.0, 0.0, 0.0]))
        # writer.write(mol)

        mol = Chem.MolFromSmiles('CO')
        rdDistGeom.EmbedMolecule(mol, 10, 1)
        conf = mol.GetConformer()
        self.assertTrue(lstEq(conf.GetAtomPosition(0), [0.69192, 0.0, 0.0]))
        self.assertTrue(lstEq(conf.GetAtomPosition(1), [-0.69192, 0.0, 0.0]))
        # writer.write(mol)

        mol = Chem.MolFromSmiles('CCC')
        rdDistGeom.EmbedMolecule(mol, 10, 1)
        conf = mol.GetConformer()
        self.assertTrue(lstEq(conf.GetAtomPosition(0), [-1.21676, -0.2989, 0.0]))
        self.assertTrue(lstEq(conf.GetAtomPosition(1), [-0.00604, 0.59337, 0.0]))
        self.assertTrue(lstEq(conf.GetAtomPosition(2), [1.22281, -0.29446, 0.0]))
        # writer.write(mol)

        mol = Chem.MolFromSmiles('O=C=O')
        rdDistGeom.EmbedMolecule(mol, 10, 1)
        conf = mol.GetConformer()

        # writer.write(mol)
        self.assertTrue(lstEq(conf.GetAtomPosition(0), [-1.2180, -0.06088, 0.0]))
        self.assertTrue(lstEq(conf.GetAtomPosition(1), [-0.00408, 0.12116, 0.0]))
        self.assertTrue(lstEq(conf.GetAtomPosition(2), [1.22207, -0.060276, 0.0]))

        mol = Chem.MolFromSmiles('C=C=C=C')
        rdDistGeom.EmbedMolecule(mol, 10, 1, useExpTorsionAnglePrefs=False,
                                 useBasicKnowledge=False)
        conf = mol.GetConformer()

        # writer.write(mol)

        d1 = computeDist(conf.GetAtomPosition(0), conf.GetAtomPosition(1))
        self.assertTrue(feq(d1, 1.31, 0.01))
        d2 = computeDist(conf.GetAtomPosition(0), conf.GetAtomPosition(2))
        self.assertTrue(feq(d2, 2.59, 0.05))
        d3 = computeDist(conf.GetAtomPosition(0), conf.GetAtomPosition(3))
        self.assertTrue(feq(d3, 3.84, 0.1))
        d4 = computeDist(conf.GetAtomPosition(1), conf.GetAtomPosition(2))
        self.assertTrue(feq(d4, 1.29, 0.01))
        d5 = computeDist(conf.GetAtomPosition(1), conf.GetAtomPosition(3))
        self.assertTrue(feq(d5, 2.54, 0.1))
        d6 = computeDist(conf.GetAtomPosition(2), conf.GetAtomPosition(3))
        self.assertTrue(feq(d6, 1.31, 0.01))

    def test2Utils(self):
        mol = Chem.MolFromSmiles('CC')
        bm = rdDistGeom.GetMoleculeBoundsMatrix(mol)
        self.assertTrue(bm[1, 0] > 0)
        self.assertTrue(bm[0, 1] > 0)
        self.assertTrue(bm[0, 1] >= bm[1, 0])
        self.assertTrue(bm[1, 0] < 1.510)
        self.assertTrue(bm[0, 1] > 1.510)

    def test3MultiConf(self):
        mol = Chem.MolFromSmiles("CC(C)(C)c(cc12)n[n]2C(=O)/C=C(N1)/COC")
        cids = rdDistGeom.EmbedMultipleConfs(mol, 10, maxAttempts=30, randomSeed=100,
                                 useExpTorsionAnglePrefs=False,
                                 useBasicKnowledge=False)
        energies = [116.330, 106.246, 109.816, 104.890,
            93.060, 140.803, 139.253, 95.820, 123.591, 108.655]
        nenergies = []
        for cid in cids:
            ff = ChemicalForceFields.UFFGetMoleculeForceField(mol, 10.0, cid)
            ee = ff.CalcEnergy()
            nenergies.append(ee)
        # print(['%.3f' % x for x in nenergies])
        # print(nenergies)
        self.assertTrue(lstEq(energies, nenergies, tol=1e-2))

    def test4OrderDependence(self):
        self.assertTrue(
          compareOrder("CC(C)(C)C(=O)NC(C1)CC(N2C)CCC12", "CN1C2CCC1CC(NC(=O)C(C)(C)C)C2"))
        # issue 230
        self.assertTrue(compareOrder("C#CC(C)(C)N(CN1)C\\N=C/1SC", "CSC1=NCN(C(C)(C)C#C)CN1"))
        # issue 232
        self.assertTrue(
          compareOrder("CC(C)(C)C(=O)NC(C1)CC(N2C)CCC12", "CN1C2CCC1CC(NC(=O)C(C)(C)C)C2"))

    def test5Issue285(self):
        m = Chem.MolFromSmiles('CNC=O')
        cs = rdDistGeom.EmbedMultipleConfs(m, 10)
        for i, ci in enumerate(cs):
            for j in range(i + 1, len(cs)):
                cj = cs[j]
                self.assertTrue(Chem.MolToMolBlock(m, confId=ci)
                                != Chem.MolToMolBlock(m, confId=cj))

    def test6RmsPruning(self):
        smiles = [
          'CC(C)CC(NC(C1[N+]CCC1)=O)C([O-])=O', 'CC(NC(CO)C(O)c1ccc([N+]([O-])=O)cc1)=O',
          'CC([N+])C(NC(C)C(N1C(C=O)CCC1)=O)=O', 'CC(NC1C(O)C=C(C([O-])=O)OC1C(O)C(O)CO)=O',
          'CCCC=C(NC(C1CC1(C)C)=O)C([O-])=O', 'OCC(O)C(O)C(Cn1c2c(cc(C)c(C)c2)nc-2c(=O)[nH]c(=O)nc12)O'
        ]

        nconfs = []
        expected = [4, 5, 5, 4, 5, 4]
        expected = [3, 3, 5, 4, 4, 4]
        for smi in smiles:
            mol = Chem.MolFromSmiles(smi)
            cids = rdDistGeom.EmbedMultipleConfs(mol, 50, maxAttempts=30, randomSeed=100,
                                                 pruneRmsThresh=1.5)
            nconfs.append(len(cids))

        d = [abs(x - y) for x, y in zip(expected, nconfs)]
        # print(nconfs)
        self.assertTrue(max(d) <= 1)

        # previous settings
        params = rdDistGeom.ETKDG()
        params.randomSeed = 100
        params.maxIterations = 30
        params.pruneRmsThresh = 1.5
        params.useSymmetryForPruning = False
        nconfs = []
        expected = [4, 5, 5, 4, 5, 4]
        for smi in smiles:
            mol = Chem.MolFromSmiles(smi)
            cids = rdDistGeom.EmbedMultipleConfs(mol, 50, params)
            nconfs.append(len(cids))

        d = [abs(x - y) for x, y in zip(expected, nconfs)]
        # print(nconfs)
        self.assertTrue(max(d) <= 1)

    def test6Chirality(self):
        # turn on chirality and we should get chiral volume that is pretty consistent and
        # positive
        tgtVol = 13.0
        smiles = "Cl[C@](C)(F)Br"
        mol = Chem.MolFromSmiles(smiles)
        cids = rdDistGeom.EmbedMultipleConfs(mol, 30, maxAttempts=30, randomSeed=100)
        self.assertTrue(len(cids) == 30)
        for cid in cids:
            conf = mol.GetConformer(cid)
            vol = computeChiralVol(
              conf.GetAtomPosition(0),
              conf.GetAtomPosition(2), conf.GetAtomPosition(3), conf.GetAtomPosition(4))
            self.assertTrue(abs(vol - tgtVol) < 1)

        # turn of chirality and now we should see both chiral forms
        smiles = "ClC(C)(F)Br"
        mol = Chem.MolFromSmiles(smiles)
        cids = rdDistGeom.EmbedMultipleConfs(mol, 30, maxAttempts=30, randomSeed=120)
        self.assertTrue(len(cids) == 30)
        nPos = 0
        nNeg = 0
        for cid in cids:
            conf = mol.GetConformer(cid)
            vol = computeChiralVol(
              conf.GetAtomPosition(0),
              conf.GetAtomPosition(2), conf.GetAtomPosition(3), conf.GetAtomPosition(4))
            self.assertTrue(abs(vol - tgtVol) < 1 or abs(vol + tgtVol) < 1)
            if vol < 0:
                nNeg += 1
            else:
                nPos += 1
        self.assertTrue(nPos > 0)
        self.assertTrue(nNeg > 0)

        tgtVol = 5.0
        for i in range(10):
            smiles = "Cl[C@H](F)Br"
            mol = Chem.MolFromSmiles(smiles)
            ci = rdDistGeom.EmbedMolecule(mol, 30, (i + 1) * 10)
            conf = mol.GetConformer(ci)
            vol = computeChiralVol(
              conf.GetAtomPosition(0),
              conf.GetAtomPosition(1), conf.GetAtomPosition(2), conf.GetAtomPosition(3))
            self.assertTrue(abs(vol - tgtVol) < 1, "%s %s" % (vol, tgtVol))

        tgtVol = 3.5
        expected = [-3.62, -3.67, -3.72, 3.91, 3.95, 3.98, 3.90, 3.94, 3.98, 3.91]
        nPos = 0
        nNeg = 0
        for i in range(30):
            smiles = "ClC(F)Br"
            mol = Chem.MolFromSmiles(smiles)
            ci = rdDistGeom.EmbedMolecule(mol, 30, (i + 1) * 10)
            conf = mol.GetConformer(ci)
            vol = computeChiralVol(
              conf.GetAtomPosition(0),
              conf.GetAtomPosition(1), conf.GetAtomPosition(2), conf.GetAtomPosition(3))
            self.assertTrue(abs(vol - tgtVol) < 1 or abs(vol + tgtVol) < 1)
            if vol < 0:
                nNeg += 1
            else:
                nPos += 1

        self.assertTrue(nPos > 0)
        self.assertTrue(nNeg > 0)

        smiles = "Cl[C@H](F)Br"
        m = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(m)
        cids = rdDistGeom.EmbedMultipleConfs(mol, 10, maxAttempts=30, randomSeed=100)
        self.assertTrue(len(cids) == 10)
        tgtVol = 10.5
        for cid in cids:
            conf = mol.GetConformer(cid)
            vol = computeChiralVol(
              conf.GetAtomPosition(0),
              conf.GetAtomPosition(2), conf.GetAtomPosition(3), conf.GetAtomPosition(4))
            self.assertTrue(abs(vol - tgtVol) < 2.)

        # let's try a little more complicated system
        expectedV1 = -2.0
        expectedV2 = -2.9

        for i in range(5):
            smi = "C1=CC=C(C=C1)[C@H](OC1=C[NH]N=C1)C(=O)[NH]C[C@H](Cl)C1=CC=NC=C1"
            mol = Chem.MolFromSmiles(smi)
            ci = rdDistGeom.EmbedMolecule(mol, randomSeed=(i + 1) * 15)
            self.assertTrue(ci >= 0)
            ff = ChemicalForceFields.UFFGetMoleculeForceField(mol, 10.0, ci)
            ff.Minimize()

            conf = mol.GetConformer(ci)
            vol1 = computeChiralVol(
              conf.GetAtomPosition(6),
              conf.GetAtomPosition(3), conf.GetAtomPosition(7), conf.GetAtomPosition(13))
            self.assertTrue(abs(vol1 - expectedV1) < 1 or abs(vol1 + expectedV1) < 1)
            if vol1 < 0:
                nNeg += 1
            else:
                nPos += 1

            vol2 = computeChiralVol(
              conf.GetAtomPosition(17),
              conf.GetAtomPosition(16), conf.GetAtomPosition(18), conf.GetAtomPosition(19))
            self.assertTrue(abs(vol2 - expectedV2) < 1 or abs(vol2 + expectedV2) < 1)

        # remove the chiral specification and we should see other chiral
        # forms of the compound
        expectedV1 = 2.0  # [-2.30, -2.31, -2.30,  2.30, -1.77]
        expectedV2 = 2.8  # [2.90,  2.89,  2.69, -2.90, -2.93]

        self.assertTrue(nPos > 0)
        self.assertTrue(nNeg > 0)
        for i in range(5):
            smi = "C1=CC=C(C=C1)C(OC1=C[NH]N=C1)C(=O)[NH]CC(Cl)C1=CC=NC=C1"
            mol = Chem.MolFromSmiles(smi)
            ci = rdDistGeom.EmbedMolecule(mol, 30, (i + 1) * 10)
            ff = ChemicalForceFields.UFFGetMoleculeForceField(mol, 10.0, ci)
            ff.Minimize()

            conf = mol.GetConformer(ci)
            vol1 = computeChiralVol(
              conf.GetAtomPosition(6),
              conf.GetAtomPosition(3), conf.GetAtomPosition(7), conf.GetAtomPosition(13))
            vol2 = computeChiralVol(
              conf.GetAtomPosition(17),
              conf.GetAtomPosition(16), conf.GetAtomPosition(18), conf.GetAtomPosition(19))
            self.assertTrue(abs(abs(vol1) - expectedV1) < 1.0)
            self.assertTrue(abs(abs(vol2) - expectedV2) < 1.0)

    def test7ConstrainedEmbedding(self):
        ofile = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'DistGeomHelpers', 'test_data',
                             'constrain1.sdf')
        suppl = Chem.SDMolSupplier(ofile)
        ref = next(suppl)
        probe = copy.deepcopy(ref)

        cMap = {}
        for i in range(5):
            cMap[i] = ref.GetConformer().GetAtomPosition(i)
        ci = rdDistGeom.EmbedMolecule(probe, coordMap=cMap, randomSeed=23)
        self.assertTrue(ci > -1)
        algMap = list(zip(range(5), range(5)))
        ssd = rdMolAlign.AlignMol(probe, ref, atomMap=algMap)
        self.assertTrue(ssd < 0.1)

    def test8MultiThreadMultiConf(self):
        if (rdBase.rdkitBuild.split('|')[2] != "MINGW"):
            ENERGY_TOLERANCE = 1.0e-6
            MSD_TOLERANCE = 1.0e-6
        else:
            ENERGY_TOLERANCE = 1.0
            MSD_TOLERANCE = 1.0e-5
        mol = Chem.AddHs(Chem.MolFromSmiles("CC(C)(C)c(cc12)n[n]2C(=O)/C=C(N1)/COC"))
        cids = rdDistGeom.EmbedMultipleConfs(mol, 200, maxAttempts=30, randomSeed=100)
        energies = []
        for cid in cids:
            ff = ChemicalForceFields.UFFGetMoleculeForceField(mol, 10.0, cid)
            ee = ff.CalcEnergy()
            energies.append(ee)

        mol2 = Chem.AddHs(Chem.MolFromSmiles("CC(C)(C)c(cc12)n[n]2C(=O)/C=C(N1)/COC"))
        cids2 = rdDistGeom.EmbedMultipleConfs(
            mol2, 200, maxAttempts=30, randomSeed=100, numThreads=4)
        self.assertTrue(lstEq(cids, cids2))
        nenergies = []
        for cid in cids2:
            ff = ChemicalForceFields.UFFGetMoleculeForceField(mol2, 10.0, cid)
            ee = ff.CalcEnergy()
            nenergies.append(ee)

        self.assertTrue(lstEq(energies, nenergies, tol=ENERGY_TOLERANCE))

        for cid in cids:
            msd = 0.0
            for i in range(mol.GetNumAtoms()):
                msd += (mol.GetConformer().GetAtomPosition(i) -
                    mol2.GetConformer().GetAtomPosition(i)).LengthSq()
            msd /= mol.GetNumAtoms()
            self.assertTrue(msd < MSD_TOLERANCE)

    def _compareConfs(self, mol, ref, molConfId, refConfId):
        self.assertEqual(mol.GetNumAtoms(), ref.GetNumAtoms())
        molConf = mol.GetConformer(molConfId)
        refConf = ref.GetConformer(refConfId)
        for i in range(mol.GetNumAtoms()):
            mp = molConf.GetAtomPosition(i)
            rp = refConf.GetAtomPosition(i)
            self.assertAlmostEqual((mp - rp).Length(), 0.0, 3)

    def test9EmbedParams(self):
        mol = Chem.AddHs(Chem.MolFromSmiles('OCCC'))
        fn = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'DistGeomHelpers', 'test_data',
                          'simple_torsion.dg.mol')
        ref = Chem.MolFromMolFile(fn, removeHs=False)
        params = rdDistGeom.EmbedParameters()
        params.randomSeed = 42
        self.assertEqual(rdDistGeom.EmbedMolecule(mol, params), 0)
        self._compareConfs(mol, ref, 0, 0)

        fn = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'DistGeomHelpers', 'test_data',
                          'simple_torsion.etdg.mol')
        ref = Chem.MolFromMolFile(fn, removeHs=False)
        params = rdDistGeom.EmbedParameters()
        params.randomSeed = 42
        params.useExpTorsionAnglePrefs = True
        self.assertEqual(rdDistGeom.EmbedMolecule(mol, params), 0)
        self._compareConfs(mol, ref, 0, 0)
        params = rdDistGeom.ETDG()
        params.randomSeed = 42
        self.assertEqual(rdDistGeom.EmbedMolecule(mol, params), 0)
        self._compareConfs(mol, ref, 0, 0)

        fn = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'DistGeomHelpers', 'test_data',
                          'simple_torsion.etkdg.mol')
        ref = Chem.MolFromMolFile(fn, removeHs=False)
        params = rdDistGeom.EmbedParameters()
        params.randomSeed = 42
        params.useExpTorsionAnglePrefs = True
        params.useBasicKnowledge = True
        self.assertEqual(rdDistGeom.EmbedMolecule(mol, params), 0)
        self._compareConfs(mol, ref, 0, 0)
        params = rdDistGeom.ETKDG()
        params.randomSeed = 42
        self.assertEqual(rdDistGeom.EmbedMolecule(mol, params), 0)
        self._compareConfs(mol, ref, 0, 0)

        fn = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'DistGeomHelpers', 'test_data',
                          'simple_torsion.kdg.mol')
        ref = Chem.MolFromMolFile(fn, removeHs=False)
        params = rdDistGeom.EmbedParameters()
        params.randomSeed = 42
        params.useBasicKnowledge = True
        self.assertEqual(rdDistGeom.EmbedMolecule(mol, params), 0)
        self._compareConfs(mol, ref, 0, 0)
        params = rdDistGeom.KDG()
        params.randomSeed = 42
        self.assertEqual(rdDistGeom.EmbedMolecule(mol, params), 0)
        self._compareConfs(mol, ref, 0, 0)

    def test10ETKDGv2(self):
        mol = Chem.AddHs(Chem.MolFromSmiles('n1cccc(C)c1ON'))
        fn = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'DistGeomHelpers', 'test_data',
                          'torsion.etkdg.v2.mol')
        ref = Chem.MolFromMolFile(fn, removeHs=False)
        params = rdDistGeom.ETKDGv2()
        params.randomSeed = 42
        self.assertEqual(rdDistGeom.EmbedMolecule(mol, params), 0)
        self._compareConfs(mol, ref, 0, 0)

    def assertDeterministicWithSeed(self, seed):
        input_mol = Chem.MolFromSmiles('CN(Cc1cnc2nc(N)nc(N)c2n1)c1ccc(C(=O)NC(CCC(=O)O)C(=O)O)cc1')

        params = AllChem.ETKDG()
        params.pruneRmsThresh = -1.0  # skip internal RMSD pruning
        if seed is not None:
            params.randomSeed = seed

        firstMol = Chem.AddHs(input_mol)
        firstIds = AllChem.EmbedMultipleConfs(firstMol, 11, params)

        secondMol = Chem.AddHs(input_mol)
        secondIds = AllChem.EmbedMultipleConfs(secondMol, 11, params)

        self.assertEqual(list(firstIds), list(secondIds))
        self.assertEqual(firstMol.GetNumConformers(), secondMol.GetNumConformers())

        nonDeterministic = False
        for confIdx in range(firstMol.GetNumConformers()):

            firstConf = firstMol.GetConformer(confIdx)
            secondConf = secondMol.GetConformer(confIdx)

            firstPositions = firstConf.GetPositions()
            secondPositions = secondConf.GetPositions()

            d = firstPositions - secondPositions
            rmsd = numpy.sqrt(numpy.sum(d * d))
            if seed >= 0:
                self.assertEqual(rmsd, 0.0)
            elif rmsd != 0.0:
                nonDeterministic = True
        if seed < 0:
            self.assertTrue(nonDeterministic)

    def testETKDGIsDeterministic(self):
        self.assertDeterministicWithSeed(-1)  # not deterministic
        self.assertDeterministicWithSeed(0)  # deterministic
        self.assertDeterministicWithSeed(1)  # deterministic
        # as large as we can go without overflowing since 11 * 195225786 should not overflow the int
        self.assertDeterministicWithSeed(195225786)
        self.assertDeterministicWithSeed(195225787)  # one higher seed will overflow though
        # another large seeds that shouldn't overflow internals and make them non-deterministic
        self.assertDeterministicWithSeed(0x1CEB00DA)

    def testGithub1763(self):
        mol = Chem.MolFromSmiles('CCCCC')
        bm1 = rdDistGeom.GetMoleculeBoundsMatrix(mol)
        bm2 = rdDistGeom.GetMoleculeBoundsMatrix(mol, doTriangleSmoothing=False)
        self.assertTrue(bm1[0, 4] < bm2[0, 4])

    def testGithub2057(self):
            # ensure that ETKDG is the default Embedder
        mol = Chem.AddHs(Chem.MolFromSmiles('OCCC'))
        fn = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'DistGeomHelpers', 'test_data',
                          'simple_torsion.etkdg.mol')
        ref = Chem.MolFromMolFile(fn, removeHs=False)
        self.assertEqual(rdDistGeom.EmbedMolecule(mol, randomSeed=42), 0)
        self._compareConfs(mol, ref, 0, 0)

    def testProvidingBoundsMatrix(self):
        m1 = Chem.MolFromSmiles("C1CCC1C")
        bm1 = rdDistGeom.GetMoleculeBoundsMatrix(m1)
        bm1[0,3] = 1.21
        bm1[3,0] = 1.20
        bm1[2,3] = 1.21
        bm1[3,2] = 1.20
        bm1[4,3] = 1.21
        bm1[3,4] = 1.20
        DG.DoTriangleSmoothing(bm1)
        ps = rdDistGeom.EmbedParameters()
        ps.useRandomCoords = True
        ps.SetBoundsMat(bm1)
        ps.randomSeed = 0xf00d
        self.assertEqual(rdDistGeom.EmbedMolecule(m1,ps),0)
        conf = m1.GetConformer()
        self.assertAlmostEqual((conf.GetAtomPosition(3)-conf.GetAtomPosition(0)).Length(),1.2,delta=0.05)
        self.assertAlmostEqual((conf.GetAtomPosition(3)-conf.GetAtomPosition(2)).Length(),1.2,delta=0.05)
        self.assertAlmostEqual((conf.GetAtomPosition(3)-conf.GetAtomPosition(4)).Length(),1.2,delta=0.05)

    def testProvidingCPCI(self):
        """
        test for a ring molecule, repeated generating a conformer with and without enforcing 
        an additional +ve interaction between a pair of non-bonded atoms (termed CPCI, 
        custom pairwise charge-like interaciton), in every iteration, applying CPCI should
        yield a conformer where this pair of atoms are further apart.
        """
        for i in range(5):
            ps = rdDistGeom.EmbedParameters()
            ps.randomSeed = i
            ps.useBasicKnowledge = True
            ps.useRandomCoords = False
            m1 = Chem.MolFromSmiles("C1CCCC1C")
            self.assertEqual(rdDistGeom.EmbedMolecule(m1,ps),0)

            m2 = Chem.MolFromSmiles("C1CCCC1C")
            ps = rdDistGeom.EmbedParameters()
            ps.randomSeed = i
            ps.useRandomCoords = False
            ps.useBasicKnowledge = True
            ps.SetCPCI({ (0,3) : 0.9 } )
            self.assertEqual(rdDistGeom.EmbedMolecule(m2,ps),0)

            conf1 = m1.GetConformer()
            conf2 = m2.GetConformer()
            self.assertTrue((conf2.GetAtomPosition(3)-conf2.GetAtomPosition(0)).Length() > (conf1.GetAtomPosition(3)-conf1.GetAtomPosition(0)).Length())

    def testScaleBoundsMatForce(self):
        """
        for pentane, set a target distance for the 1-5 distance, and generate conformers with changing weights for (all) the atom pair distance restraints,
        the conformer with the stronger weight for the atom pairs will always have a 1-5 distance closer to the target value than that with the weaker weight.
        """
        target = 4
        for i in range(5):
            ps = rdDistGeom.EmbedParameters()
            ps.randomSeed = i
            ps.useBasicKnowledge = True
            ps.useRandomCoords = False
            m1 = Chem.MolFromSmiles("CCCCC")
            bm1 = rdDistGeom.GetMoleculeBoundsMatrix(m1)
            bm1[0,4] = target
            bm1[4,0] = target
            DG.DoTriangleSmoothing(bm1)
            ps.boundsMatForceScaling = 0.1
            ps.SetBoundsMat(bm1)
            self.assertEqual(rdDistGeom.EmbedMolecule(m1,ps),0)

            m2 = Chem.MolFromSmiles("CCCCC")
            ps = rdDistGeom.EmbedParameters()
            ps.randomSeed = i
            ps.useBasicKnowledge = True
            ps.useRandomCoords = False
            ps.boundsMatForceScaling = 10
            ps.SetBoundsMat(bm1)
            self.assertEqual(rdDistGeom.EmbedMolecule(m2,ps),0)

            conf1 = m1.GetConformer()
            conf2 = m2.GetConformer()
            self.assertTrue(abs((conf2.GetAtomPosition(4)-conf2.GetAtomPosition(0)).Length() - target) < abs((conf1.GetAtomPosition(4)-conf1.GetAtomPosition(0)).Length() - target))

        
    def testETKDGv3amide(self):
        """
        test for a macrocycle molecule, ETKDGv3 samples trans amide
        """
        def get_atom_mapping(mol, smirks = "[O:1]=[C:2]@;-[NX3:3]-[H:4]"):
            qmol = Chem.MolFromSmarts(smirks)
            ind_map = {}
            for atom in qmol.GetAtoms():
                map_num = atom.GetAtomMapNum()
                if map_num:
                    ind_map[map_num - 1] = atom.GetIdx()
            map_list = [ind_map[x] for x in sorted(ind_map)]
            matches = list()
            for match in mol.GetSubstructMatches(qmol, uniquify = False) :
                mas = [match[x] for x in map_list]
                matches.append(tuple(mas))
            return matches

        smiles = "C1CCC(=O)NCCCCCC(=O)NC1"
        smiles_mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(smiles_mol)
        params = AllChem.ETKDGv3()
        params.seed = 42
        AllChem.EmbedMolecule(mol, params)
        conf = mol.GetConformer(0)
        for torsion in get_atom_mapping(mol):
            a1,a2,a3,a4 = [conf.GetAtomPosition(i) for i in torsion]
            self.assertAlmostEqual(abs(ComputeSignedDihedralAngle(a1,a2,a3,a4)), 3.14, delta = 0.1)


if __name__ == '__main__':
    unittest.main()
