import unittest
from rdkit import Chem
from rdkit.Chem import rdFMCS


class BondMatchOrderMatrix:
    def __init__(self, ignoreAromatization):
        self.MatchMatrix = [[False]*(Chem.BondType.ZERO + 1)
                            for i in range(Chem.BondType.ZERO + 1)]
        for i in range(Chem.BondType.ZERO + 1):
            # fill cells of the same and unspecified type
            self.MatchMatrix[i][i] = True
            self.MatchMatrix[Chem.BondType.UNSPECIFIED][i] = \
                self.MatchMatrix[i][Chem.BondType.UNSPECIFIED] = True
            self.MatchMatrix[Chem.BondType.ZERO][i] = \
                self.MatchMatrix[i][Chem.BondType.ZERO] = True
        if ignoreAromatization:
            self.MatchMatrix[Chem.BondType.SINGLE][Chem.BondType.AROMATIC] = \
                self.MatchMatrix[Chem.BondType.AROMATIC][Chem.BondType.SINGLE] = True
            self.MatchMatrix[Chem.BondType.DOUBLE][Chem.BondType.AROMATIC] = \
                self.MatchMatrix[Chem.BondType.AROMATIC][Chem.BondType.DOUBLE] = True
        self.MatchMatrix[Chem.BondType.SINGLE][Chem.BondType.ONEANDAHALF] = \
            self.MatchMatrix[Chem.BondType.ONEANDAHALF][Chem.BondType.SINGLE] = True
        self.MatchMatrix[Chem.BondType.DOUBLE][Chem.BondType.TWOANDAHALF] = \
            self.MatchMatrix[Chem.BondType.TWOANDAHALF][Chem.BondType.DOUBLE] = True
        self.MatchMatrix[Chem.BondType.TRIPLE][Chem.BondType.THREEANDAHALF] = \
            self.MatchMatrix[Chem.BondType.THREEANDAHALF][Chem.BondType.TRIPLE] = True
        self.MatchMatrix[Chem.BondType.QUADRUPLE][Chem.BondType.FOURANDAHALF] = \
            self.MatchMatrix[Chem.BondType.FOURANDAHALF][Chem.BondType.QUADRUPLE] = True
        self.MatchMatrix[Chem.BondType.QUINTUPLE][Chem.BondType.FIVEANDAHALF] = \
            self.MatchMatrix[Chem.BondType.FIVEANDAHALF][Chem.BondType.QUINTUPLE] = True
    def isEqual(self, i, j):
        return self.MatchMatrix[i][j]

class CompareAny(rdFMCS.MCSAtomCompare):
    def compare(self, p, mol1, atom1, mol2, atom2):
        if (p.MatchChiralTag and not self.CheckAtomChirality(p, mol1, atom1, mol2, atom2)):
            return False
        if (p.MatchFormalCharge and not self.CheckAtomCharge(p, mol1, atom1, mol2, atom2)):
            return False
        if (p.RingMatchesRingOnly):
            return self.CheckAtomRingMatch(p, mol1, atom1, mol2, atom2)
        return True

class CompareAnyHeavyAtom(CompareAny):
    def compare(self, p, mol1, atom1, mol2, atom2):
        a1 = mol1.GetAtomWithIdx(atom1)
        a2 = mol2.GetAtomWithIdx(atom2)
        # Any atom, including H, matches another atom of the same type,  according to
        # the other flags
        if (a1.GetAtomicNum() == a2.GetAtomicNum() or
            (a1.GetAtomicNum() > 1 and a2.GetAtomicNum() > 1)):
            return CompareAny.compare(self, p, mol1, atom1, mol2, atom2)
        return False

class CompareElements(rdFMCS.MCSAtomCompare):
    def compare(self, p, mol1, atom1, mol2, atom2):
        a1 = mol1.GetAtomWithIdx(atom1)
        a2 = mol2.GetAtomWithIdx(atom2)
        if (a1.GetAtomicNum() != a2.GetAtomicNum()):
            return False
        if (p.MatchValences and a1.GetTotalValence() != a2.GetTotalValence()):
            return False
        if (p.MatchChiralTag and not self.CheckAtomChirality(p, mol1, atom1, mol2, atom2)):
            return False
        if (p.MatchFormalCharge and not self.CheckAtomCharge(p, mol1, atom1, mol2, atom2)):
            return False
        if p.RingMatchesRingOnly:
            return self.CheckAtomRingMatch(p, mol1, atom1, mol2, atom2)
        return True

class CompareIsotopes(rdFMCS.MCSAtomCompare):
    def compare(self, p, mol1, atom1, mol2, atom2):
        a1 = mol1.GetAtomWithIdx(atom1)
        a2 = mol2.GetAtomWithIdx(atom2)
        if (a1.GetIsotope() != a2.GetIsotope()):
            return False
        if (p.MatchChiralTag and not self.CheckAtomChirality(p, mol1, atom1, mol2, atom2)):
            return False
        if (p.MatchFormalCharge and not self.CheckAtomCharge(p, mol1, atom1, mol2, atom2)):
            return False
        if p.RingMatchesRingOnly:
            return self.CheckAtomRingMatch(p, mol1, atom1, mol2, atom2)
        return True

class CompareOrder(rdFMCS.MCSBondCompare):
    match = BondMatchOrderMatrix(True)  # ignore Aromatization
    def compare(self, p, mol1, bond1, mol2, bond2):
        b1 = mol1.GetBondWithIdx(bond1)
        b2 = mol2.GetBondWithIdx(bond2)
        t1 = b1.GetBondType()
        t2 = b2.GetBondType()
        if self.match.isEqual(t1, t2):
            if (p.MatchStereo and not self.CheckBondStereo(p, mol1, bond1, mol2, bond2)):
                return False
            if p.RingMatchesRingOnly:
                return self.CheckBondRingMatch(p, mol1, bond1, mol2, bond2)
            return True
        return False

class AtomCompareCompareIsInt(rdFMCS.MCSAtomCompare):
    compare = 1

class AtomCompareNoCompare(rdFMCS.MCSAtomCompare):
    pass

class AtomCompareUserData(rdFMCS.MCSAtomCompare):
    def __init__(self):
        super().__init__()
        self._matchAnyHet = False
    def setMatchAnyHet(self, v):
        self._matchAnyHet = v
    def compare(self, p, mol1, atom1, mol2, atom2):
        a1 = mol1.GetAtomWithIdx(atom1)
        a2 = mol2.GetAtomWithIdx(atom2)
        if (a1.GetAtomicNum() != a2.GetAtomicNum() and
            ((not self._matchAnyHet) or
            a1.GetAtomicNum() == 6 or
            a2.GetAtomicNum() == 6)):
            return False
        if (p.MatchValences and a1.GetTotalValence() != a2.GetTotalValence()):
            return False
        if (p.MatchChiralTag and not self.CheckAtomChirality(p, mol1, atom1, mol2, atom2)):
            return False
        if (p.MatchFormalCharge and not self.CheckAtomCharge(p, mol1, atom1, mol2, atom2)):
            return False
        if p.RingMatchesRingOnly:
            return self.CheckAtomRingMatch(p, mol1, atom1, mol2, atom2)
        return True

class BondCompareCompareIsInt(rdFMCS.MCSBondCompare):
    compare = 1

class BondCompareNoCompare(rdFMCS.MCSBondCompare):
    pass

class BondCompareUserData(rdFMCS.MCSBondCompare):
    def __init__(self):
        super().__init__()
        self.match = None
    def setIgnoreAromatization(self, v):
        self.match = BondMatchOrderMatrix(v)
    def compare(self, p, mol1, bond1, mol2, bond2):
        b1 = mol1.GetBondWithIdx(bond1)
        b2 = mol2.GetBondWithIdx(bond2)
        t1 = b1.GetBondType()
        t2 = b2.GetBondType()
        if self.match.isEqual(t1, t2):
            if (p.MatchStereo and not self.CheckBondStereo(p, mol1, bond1, mol2, bond2)):
                return False
            if p.RingMatchesRingOnly:
                return self.CheckBondRingMatch(p, mol1, bond1, mol2, bond2)
            return True
        return False

class ProgressCallbackCallbackIsInt(rdFMCS.MCSProgress):
    callback = 1

class ProgressCallbackNoCallback(rdFMCS.MCSProgress):
    pass

class ProgressCallback(rdFMCS.MCSProgress):
    def __init__(self, parent):
        super().__init__()
        self.parent = parent
        self.callCount = 0
    def callback(self, stat, params):
        self.callCount += 1
        self.parent.assertTrue(isinstance(stat, rdFMCS.MCSProgressData))
        self.parent.assertTrue(hasattr(stat, "numAtoms"))
        self.parent.assertTrue(isinstance(stat.numAtoms, int))
        self.parent.assertTrue(hasattr(stat, "numBonds"))
        self.parent.assertTrue(isinstance(stat.numBonds, int))
        self.parent.assertTrue(hasattr(stat, "seedProcessed"))
        self.parent.assertTrue(isinstance(stat.seedProcessed, int))
        self.parent.assertTrue(isinstance(params, rdFMCS.MCSParameters))
        self.parent.assertTrue(isinstance(params.AtomTyper, rdFMCS.MCSAtomCompare))
        self.parent.assertTrue(isinstance(params.BondTyper, rdFMCS.BondCompare))
        self.parent.assertEqual(params.ProgressCallback, self)
        return (self.callCount < 3)

class Common:
    @staticmethod
    def getParams(**kwargs):
        have_kw = False
        params = rdFMCS.MCSParameters()
        for kw in ("AtomTyper", "BondTyper"):
            try:
                v = kwargs[kw]
            except KeyError:
                pass
            else:
                have_kw = True
                setattr(params, kw, v())
        return params

    @staticmethod
    def test1(self, **kwargs):
        smis = (
          "Cc1nc(CN(C(C)c2ncccc2)CCCCN)ccc1 CHEMBL1682991",  # -- QUERY
          "Cc1ccc(CN(C(C)c2ccccn2)CCCCN)nc1 CHEMBL1682990",
          "Cc1cccnc1CN(C(C)c1ccccn1)CCCCN CHEMBL1682998",
          "CC(N(CCCCN)Cc1c(N)cccn1)c1ccccn1 CHEMBL1682987",
          "Cc1cc(C)c(CN(C(C)c2ccccn2)CCCCN)nc1 CHEMBL1682992",
          "Cc1cc(C(C)N(CCCCN)Cc2c(C)cccn2)ncc1 CHEMBL1682993",
          "Cc1nc(C(C)N(CCCCN)Cc2nc3c([nH]2)cccc3)ccc1 CHEMBL1682878",
          "CC(c1ncccc1)N(CCCCN)Cc1nc2c([nH]1)cccc2 CHEMBL1682867",
          "CC(N(CCCCN)Cc1c(C(C)(C)C)cccn1)c1ccccn1 CHEMBL1682989",
          "CC(N(CCCCN)Cc1c(C(F)(F)F)cccn1)c1ccccn1 CHEMBL1682988",
        )
        ms = [Chem.MolFromSmiles(x.split()[0]) for x in smis]
        qm = ms[0]
        ms = ms[1:]
        if kwargs:
            params = Common.getParams(**kwargs)
            mcs = rdFMCS.FindMCS(ms, params)
        else:
            mcs = rdFMCS.FindMCS(ms)
        self.assertEqual(mcs.numBonds, 21)
        self.assertEqual(mcs.numAtoms, 21)
        self.assertEqual(
          mcs.smartsString,
          '[#6](:[#6]:[#6]):[#6]:[#7]:[#6]-[#6]-[#7](-[#6](-[#6])-[#6]1:[#6]:[#6]:[#6]:[#6]:[#7]:1)-[#6]-[#6]-[#6]-[#6]-[#7]'
        )
        qm = Chem.MolFromSmarts(mcs.smartsString)
        self.assertTrue(qm is not None)
        for m in ms:
            self.assertTrue(m.HasSubstructMatch(qm))

        if kwargs:
            params = Common.getParams(**kwargs)
            params.Threshold = 0.8
            mcs = rdFMCS.FindMCS(ms, params)
        else:
            mcs = rdFMCS.FindMCS(ms, threshold=0.8)
        self.assertEqual(mcs.numBonds, 21)
        self.assertEqual(mcs.numAtoms, 21)
        self.assertEqual(
          mcs.smartsString,
          '[#6](:[#6]:[#6]):[#6]:[#7]:[#6]-[#6]-[#7](-[#6](-[#6])-[#6]1:[#6]:[#6]:[#6]:[#6]:[#7]:1)-[#6]-[#6]-[#6]-[#6]-[#7]'
        )
        qm = Chem.MolFromSmarts(mcs.smartsString)
        self.assertTrue(qm is not None)
        for m in ms:
            self.assertTrue(m.HasSubstructMatch(qm))

    def test2(self, **kwargs):
        smis = (
          "CHEMBL122452 CN(CCCN(C)CCc1ccccc1)CCOC(c1ccccc1)c1ccccc1",
          "CHEMBL123252 CN(CCCc1ccccc1)CCCN(C)CCOC(c1ccccc1)c1ccccc1",
          "CHEMBL121611 Fc1ccc(C(OCCNCCCNCCc2ccccc2)c2ccc(F)cc2)cc1",
          "CHEMBL121050 O=C(Cc1ccccc1)NCCCCNCCOC(c1ccc(F)cc1)c1ccc(F)cc1",
          "CHEMBL333667 O=C(Cc1ccccc1)NCCNCCOC(c1ccc(F)cc1)c1ccc(F)cc1",
          "CHEMBL121486 O=C(Cc1ccc(Br)cc1)NC=CNCCOC(c1ccc(F)cc1)c1ccc(F)cc1",
          "CHEMBL123830 O=C(Cc1ccc(F)cc1)NCCNCCOC(c1ccc(F)cc1)c1ccc(F)cc1",
          "CHEMBL420900 O=C(Cc1ccccc1)NCCCNCCOC(c1ccc(F)cc1)c1ccc(F)cc1",
          "CHEMBL121460 CN(CCOC(c1ccc(F)cc1)c1ccc(F)cc1)CCN(C)CCOC(c1ccc(F)cc1)c1ccc(F)cc1",
          "CHEMBL120901 COC(=O)C1C2CCC(CC1C(=O)Oc1ccccc1)N2C",
          "CHEMBL122859 O=C1CN(CCc2ccccc2)CCN1CCOC(c1ccc(F)cc1)c1ccc(F)cc1",
          "CHEMBL121027 CN(CCOC(c1ccccc1)c1ccccc1)CCN(C)CCc1ccc(F)cc1",
        )

        ms = [Chem.MolFromSmiles(x.split()[1]) for x in smis]
        qm = ms[0]
        ms = ms[1:]
        if kwargs:
            params = Common.getParams(**kwargs)
            mcs = rdFMCS.FindMCS(ms, params)
        else:
            mcs = rdFMCS.FindMCS(ms)
        self.assertEqual(mcs.numBonds, 9)
        self.assertEqual(mcs.numAtoms, 10)
        qm = Chem.MolFromSmarts(mcs.smartsString)
        self.assertTrue(qm is not None)
        for m in ms:
            self.assertTrue(m.HasSubstructMatch(qm))
        # smarts too hard to canonicalize this
        # self.assertEqual(mcs.smartsString,'[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6](-[#6]-[#8]-[#6]:,-[#6])-,:[#6]')

        if kwargs:
            params = Common.getParams(**kwargs)
            params.Threshold = 0.8
            mcs = rdFMCS.FindMCS(ms, params)
        else:
            mcs = rdFMCS.FindMCS(ms, threshold=0.8)
        self.assertEqual(mcs.numBonds, 20)
        self.assertEqual(mcs.numAtoms, 19)
        qm = Chem.MolFromSmarts(mcs.smartsString)
        self.assertTrue(qm is not None)
        nHits = 0
        for m in ms:
            if m.HasSubstructMatch(qm):
                nHits += 1
        self.assertTrue(nHits >= int(0.8 * len(smis)))
        # smarts too hard to canonicalize this
        # self.assertEqual(mcs.smartsString,'[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#6](-[#8]-[#6]-[#6]-[#7]-[#6]-[#6])-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2')

    def test3IsotopeMatch(self, **kwargs):
        smis = (
          "CC[14NH2]",
          "CC[14CH3]",
        )

        ms = [Chem.MolFromSmiles(x) for x in smis]
        if kwargs:
            params = Common.getParams(**kwargs)
            mcs = rdFMCS.FindMCS(ms, params)
        else:
            mcs = rdFMCS.FindMCS(ms)
        self.assertEqual(mcs.numBonds, 1)
        self.assertEqual(mcs.numAtoms, 2)
        qm = Chem.MolFromSmarts(mcs.smartsString)

        if kwargs:
            params = Common.getParams(**kwargs)
            params.AtomTyper = CompareIsotopes()
            params.AtomCompareParameters.MatchIsotope = True
            mcs = rdFMCS.FindMCS(ms, params)
        else:
            mcs = rdFMCS.FindMCS(ms, atomCompare=rdFMCS.AtomCompare.CompareIsotopes)
        self.assertEqual(mcs.numBonds, 2)
        self.assertEqual(mcs.numAtoms, 3)
        qm = Chem.MolFromSmarts(mcs.smartsString)

        self.assertTrue(Chem.MolFromSmiles('CC[14CH3]').HasSubstructMatch(qm))
        self.assertFalse(Chem.MolFromSmiles('CC[13CH3]').HasSubstructMatch(qm))
        self.assertTrue(Chem.MolFromSmiles('OO[14CH3]').HasSubstructMatch(qm))
        self.assertFalse(Chem.MolFromSmiles('O[13CH2][14CH3]').HasSubstructMatch(qm))

    def test4RingMatches(self, **kwargs):
        smis = ['CCCCC', 'CCC1CCCCC1']
        ms = [Chem.MolFromSmiles(x) for x in smis]
        if kwargs:
            params = Common.getParams(**kwargs)
            mcs = rdFMCS.FindMCS(ms, params)
        else:
            mcs = rdFMCS.FindMCS(ms)
        self.assertEqual(mcs.numBonds, 4)
        self.assertEqual(mcs.numAtoms, 5)
        self.assertEqual(mcs.smartsString, '[#6]-[#6]-[#6]-[#6]-[#6]')

        if kwargs:
            params = Common.getParams(**kwargs)
            params.BondCompareParameters.CompleteRingsOnly = True
            mcs = rdFMCS.FindMCS(ms, params)
        else:
            mcs = rdFMCS.FindMCS(ms, completeRingsOnly=True)
        self.assertEqual(mcs.numBonds, 2)
        self.assertEqual(mcs.numAtoms, 3)
        self.assertEqual(mcs.smartsString, '[#6]-&!@[#6]-&!@[#6]')

        if kwargs:
            params = Common.getParams(**kwargs)
            params.BondCompareParameters.CompleteRingsOnly = True
            params.BondCompareParameters.MatchFusedRings = True
            params.BondCompareParameters.MatchFusedRingsStrict = False
            mcs = rdFMCS.FindMCS(ms, params)
        else:
            mcs = rdFMCS.FindMCS(ms, completeRingsOnly=True,
                                 ringCompare=rdFMCS.RingCompare.PermissiveRingFusion)
        self.assertEqual(mcs.numBonds, 2)
        self.assertEqual(mcs.numAtoms, 3)
        self.assertEqual(mcs.smartsString, '[#6]-&!@[#6]-&!@[#6]')

        if kwargs:
            params = Common.getParams(**kwargs)
            params.BondCompareParameters.CompleteRingsOnly = True
            params.BondCompareParameters.MatchFusedRings = True
            params.BondCompareParameters.MatchFusedRingsStrict = True
            mcs = rdFMCS.FindMCS(ms, params)
        else:
            mcs = rdFMCS.FindMCS(ms, completeRingsOnly=True,
                                 ringCompare=rdFMCS.RingCompare.StrictRingFusion)
        self.assertEqual(mcs.numBonds, 2)
        self.assertEqual(mcs.numAtoms, 3)
        self.assertEqual(mcs.smartsString, '[#6]-&!@[#6]-&!@[#6]')

        if kwargs:
            params = Common.getParams(**kwargs)
            params.AtomCompareParameters.RingMatchesRingOnly = True
            params.BondCompareParameters.RingMatchesRingOnly = True
            mcs = rdFMCS.FindMCS(ms, params)
        else:
            mcs = rdFMCS.FindMCS(ms, ringMatchesRingOnly=True)
        self.assertEqual(mcs.numBonds, 1)
        self.assertEqual(mcs.numAtoms, 2)
        self.assertEqual(mcs.smartsString, '[#6&!R]-&!@[#6&!R]')

        smis = ['CC1CCC1', 'CCC1CCCCC1']
        ms = [Chem.MolFromSmiles(x) for x in smis]
        if kwargs:
            params = Common.getParams(**kwargs)
            mcs = rdFMCS.FindMCS(ms, params)
        else:
            mcs = rdFMCS.FindMCS(ms)
        self.assertEqual(mcs.numBonds, 4)
        self.assertEqual(mcs.numAtoms, 5)
        self.assertEqual(mcs.smartsString, '[#6]-[#6](-[#6]-[#6])-[#6]')

        if kwargs:
            params = Common.getParams(**kwargs)
            params.BondCompareParameters.CompleteRingsOnly = True
            mcs = rdFMCS.FindMCS(ms, params)
        else:
            mcs = rdFMCS.FindMCS(ms, completeRingsOnly=True)
        self.assertEqual(mcs.numBonds, 1)
        self.assertEqual(mcs.numAtoms, 2)
        self.assertEqual(mcs.smartsString, '[#6]-&!@[#6]')

        if kwargs:
            params = Common.getParams(**kwargs)
            params.AtomCompareParameters.RingMatchesRingOnly = True
            params.BondCompareParameters.CompleteRingsOnly = True
            params.BondCompareParameters.RingMatchesRingOnly = True
            mcs = rdFMCS.FindMCS(ms, params)
        else:
            mcs = rdFMCS.FindMCS(ms, ringMatchesRingOnly=True, completeRingsOnly=True)
        self.assertEqual(mcs.numBonds, 1)
        self.assertEqual(mcs.numAtoms, 2)
        self.assertEqual(mcs.smartsString, '[#6&!R]-&!@[#6&R]')

        if kwargs:
            params = Common.getParams(**kwargs)
            params.AtomCompareParameters.RingMatchesRingOnly = True
            params.BondCompareParameters.CompleteRingsOnly = True
            params.BondCompareParameters.RingMatchesRingOnly = True
            params.BondCompareParameters.MatchFusedRings = True
            params.BondCompareParameters.MatchFusedRingsStrict = False
            mcs = rdFMCS.FindMCS(ms, params)
        else:
            mcs = rdFMCS.FindMCS(ms, ringMatchesRingOnly=True, completeRingsOnly=True,
                                 ringCompare=rdFMCS.RingCompare.PermissiveRingFusion)
        self.assertEqual(mcs.numBonds, 1)
        self.assertEqual(mcs.numAtoms, 2)
        self.assertEqual(mcs.smartsString, '[#6&!R]-&!@[#6&R]')

        if kwargs:
            params = Common.getParams(**kwargs)
            params.AtomCompareParameters.RingMatchesRingOnly = True
            params.BondCompareParameters.CompleteRingsOnly = True
            params.BondCompareParameters.RingMatchesRingOnly = True
            params.BondCompareParameters.MatchFusedRings = True
            params.BondCompareParameters.MatchFusedRingsStrict = True
            mcs = rdFMCS.FindMCS(ms, params)
        else:
            mcs = rdFMCS.FindMCS(ms, ringMatchesRingOnly=True, completeRingsOnly=True,
                                 ringCompare=rdFMCS.RingCompare.StrictRingFusion)
        self.assertEqual(mcs.numBonds, 1)
        self.assertEqual(mcs.numAtoms, 2)
        self.assertEqual(mcs.smartsString, '[#6&!R]-&!@[#6&R]')

        if kwargs:
            params = Common.getParams(**kwargs)
            params.AtomCompareParameters.RingMatchesRingOnly = True
            params.BondCompareParameters.RingMatchesRingOnly = True
            mcs = rdFMCS.FindMCS(ms, params)
        else:
            mcs = rdFMCS.FindMCS(ms, ringMatchesRingOnly=True)
        self.assertEqual(mcs.numBonds, 4)
        self.assertEqual(mcs.numAtoms, 5)
        self.assertEqual(mcs.smartsString, '[#6&!R]-&!@[#6&R](-&@[#6&R]-&@[#6&R])-&@[#6&R]')

    def test5AnyMatch(self, **kwargs):
        smis = ('c1ccccc1C', 'c1ccccc1O', 'c1ccccc1Cl')
        ms = [Chem.MolFromSmiles(x) for x in smis]
        if kwargs:
            params = Common.getParams(**kwargs)
            params.AtomTyper = CompareAny()
            mcs = rdFMCS.FindMCS(ms, params)
        else:
            mcs = rdFMCS.FindMCS(ms, atomCompare=rdFMCS.AtomCompare.CompareAny)
        self.assertEqual(mcs.numBonds, 7)
        self.assertEqual(mcs.numAtoms, 7)
        qm = Chem.MolFromSmarts(mcs.smartsString)

        for m in ms:
            self.assertTrue(m.HasSubstructMatch(qm))

        smis = ('c1cccnc1C', 'c1cnncc1O', 'c1cccnc1Cl')
        ms = [Chem.MolFromSmiles(x) for x in smis]
        if kwargs:
            params = Common.getParams(**kwargs)
            params.AtomTyper = CompareAny()
            mcs = rdFMCS.FindMCS(ms, params)
        else:
            mcs = rdFMCS.FindMCS(ms, atomCompare=rdFMCS.AtomCompare.CompareAny)
        self.assertEqual(mcs.numBonds, 7)
        self.assertEqual(mcs.numAtoms, 7)
        qm = Chem.MolFromSmarts(mcs.smartsString)

        for m in ms:
            self.assertTrue(m.HasSubstructMatch(qm))

    def testAtomCompareAnyHeavyAtom(self, **kwargs):
        # H matches H, O matches C
        smis = ('[H]c1ccccc1C', '[H]c1ccccc1O')
        ms = [Chem.MolFromSmiles(x, sanitize=False) for x in smis]
        if kwargs:
            params = Common.getParams(**kwargs)
            params.AtomTyper = CompareAnyHeavyAtom()
            mcs = rdFMCS.FindMCS(ms, params)
        else:
            mcs = rdFMCS.FindMCS(ms, atomCompare=rdFMCS.AtomCompare.CompareAnyHeavyAtom)
        self.assertEqual(mcs.numBonds, 8)
        self.assertEqual(mcs.numAtoms, 8)
        qm = Chem.MolFromSmarts(mcs.smartsString)

        for m in ms:
            self.assertTrue(m.HasSubstructMatch(qm))

    def testAtomCompareAnyHeavyAtom1(self, **kwargs):
        # O matches C, H does not match O
        smis = ('[H]c1ccccc1C', 'Oc1ccccc1O')
        ms = [Chem.MolFromSmiles(x, sanitize=False) for x in smis]
        if kwargs:
            params = Common.getParams(**kwargs)
            params.AtomTyper = CompareAnyHeavyAtom()
            mcs = rdFMCS.FindMCS(ms, params)
        else:
            mcs = rdFMCS.FindMCS(ms, atomCompare=rdFMCS.AtomCompare.CompareAnyHeavyAtom)
        self.assertEqual(mcs.numBonds, 7)
        self.assertEqual(mcs.numAtoms, 7)
        qm = Chem.MolFromSmarts(mcs.smartsString)

        for m in ms:
            self.assertTrue(m.HasSubstructMatch(qm))

    def test6MatchValences(self, **kwargs):
        ms = (Chem.MolFromSmiles('NC1OC1'), Chem.MolFromSmiles('C1OC1[N+](=O)[O-]'))
        if kwargs:
            params = Common.getParams(**kwargs)
            mcs = rdFMCS.FindMCS(ms, params)
        else:
            mcs = rdFMCS.FindMCS(ms)
        self.assertEqual(mcs.numBonds, 4)
        self.assertEqual(mcs.numAtoms, 4)
        if kwargs:
            params = Common.getParams(**kwargs)
            params.AtomCompareParameters.MatchValences = True
            mcs = rdFMCS.FindMCS(ms, params)
        else:
            mcs = rdFMCS.FindMCS(ms, matchValences=True)
        self.assertEqual(mcs.numBonds, 3)
        self.assertEqual(mcs.numAtoms, 3)

    def test7Seed(self, **kwargs):
        smis = ['C1CCC1CC1CC1', 'C1CCC1OC1CC1', 'C1CCC1NC1CC1', 'C1CCC1SC1CC1']
        ms = [Chem.MolFromSmiles(x) for x in smis]
        if kwargs:
            params = Common.getParams(**kwargs)
            r = rdFMCS.FindMCS(ms, params)
        else:
            r = rdFMCS.FindMCS(ms)
        self.assertEqual(r.smartsString, "[#6]1-[#6]-[#6]-[#6]-1")
        if kwargs:
            params = Common.getParams(**kwargs)
            params.InitialSeed = 'C1CC1'
            r = rdFMCS.FindMCS(ms, params)
        else:
            r = rdFMCS.FindMCS(ms, seedSmarts='C1CC1')
        self.assertEqual(r.smartsString, "[#6]1-[#6]-[#6]-1")
        if kwargs:
            params = Common.getParams(**kwargs)
            params.InitialSeed = 'C1OC1'
            r = rdFMCS.FindMCS(ms, params)
        else:
            r = rdFMCS.FindMCS(ms, seedSmarts='C1OC1')
        self.assertEqual(r.smartsString, "")

    def test8MatchParams(self, **kwargs):
        smis = ("CCC1NC1", "CCC1N(C)C1", "CCC1OC1")
        ms = [Chem.MolFromSmiles(x) for x in smis]

        if kwargs:
            ps = Common.getParams(**kwargs)
            mcs = rdFMCS.FindMCS(ms, ps)
        else:
            mcs = rdFMCS.FindMCS(ms)
        self.assertEqual(mcs.numAtoms, 4)

        ps = Common.getParams(**kwargs)
        ps.BondCompareParameters.CompleteRingsOnly = True
        mcs = rdFMCS.FindMCS(ms, ps)
        self.assertEqual(mcs.numAtoms, 3)

        ps = Common.getParams(**kwargs)
        ps.BondCompareParameters.CompleteRingsOnly = True;
        ps.BondCompareParameters.MatchFusedRings = True
        mcs = rdFMCS.FindMCS(ms, ps)
        self.assertEqual(mcs.numAtoms, 3)

        ps = Common.getParams(**kwargs)
        ps.BondCompareParameters.CompleteRingsOnly = True;
        ps.BondCompareParameters.MatchFusedRingsStrict = True
        mcs = rdFMCS.FindMCS(ms, ps)
        self.assertEqual(mcs.numAtoms, 3)

        ps = Common.getParams(**kwargs)
        if kwargs:
            ps.SetAtomTyper(CompareAny())
        else:
            ps.SetAtomTyper(rdFMCS.AtomCompare.CompareAny)
        mcs = rdFMCS.FindMCS(ms, ps)
        self.assertEqual(mcs.numAtoms, 5)

    def test9MatchCharge(self, **kwargs):
        smis = ("CCNC", "CCN(C)C", "CC[N+](C)C")
        ms = [Chem.MolFromSmiles(x) for x in smis]

        if kwargs:
            ps = Common.getParams(**kwargs)
            mcs = rdFMCS.FindMCS(ms, ps)
        else:
            mcs = rdFMCS.FindMCS(ms)

        self.assertEqual(mcs.numAtoms, 4)

        ps = Common.getParams(**kwargs)
        ps.AtomCompareParameters.MatchFormalCharge = True
        mcs = rdFMCS.FindMCS(ms, ps)
        self.assertEqual(mcs.numAtoms, 2)

    def test10MatchChargeAndParams(self, **kwargs):
        smis = ("CCNC", "CCN(C)C", "CC[N+](C)C", "CC[C+](C)C")
        ms = [Chem.MolFromSmiles(x) for x in smis]

        if kwargs:
            ps = Common.getParams(**kwargs)
            mcs = rdFMCS.FindMCS(ms, ps)
        else:
            mcs = rdFMCS.FindMCS(ms)
        self.assertEqual(mcs.numAtoms, 2)

        ps = Common.getParams(**kwargs)
        if kwargs:
            ps.AtomTyper = CompareAny()
        else:
            ps.AtomTyper = rdFMCS.AtomCompare.CompareAny
        mcs = rdFMCS.FindMCS(ms, ps)
        self.assertEqual(mcs.numAtoms, 4)

        ps = Common.getParams(**kwargs)
        ps.AtomCompareParameters.MatchFormalCharge = True
        mcs = rdFMCS.FindMCS(ms, ps)
        self.assertEqual(mcs.numAtoms, 2)

    def test11Github2034(self, **kwargs):
        smis = ("C1CC1N2CC2", "C1CC1N")
        ms = [Chem.MolFromSmiles(x) for x in smis]

        if kwargs:
            ps = Common.getParams(**kwargs)
            mcs = rdFMCS.FindMCS(ms, ps)
        else:
            mcs = rdFMCS.FindMCS(ms)
        self.assertEqual(mcs.numAtoms, 4)
        self.assertEqual(mcs.numBonds, 4)

        if kwargs:
            ps = Common.getParams(**kwargs)
            ps.AtomCompareParameters.RingMatchesRingOnly = True
            ps.BondCompareParameters.RingMatchesRingOnly = True
            mcs = rdFMCS.FindMCS(ms, ps)
        else:
            mcs = rdFMCS.FindMCS(ms, ringMatchesRingOnly=True)
        self.assertEqual(mcs.numAtoms, 3)
        self.assertEqual(mcs.numBonds, 3)

        ps = Common.getParams(**kwargs)
        ps.AtomCompareParameters.RingMatchesRingOnly = True
        mcs = rdFMCS.FindMCS(ms, ps)
        self.assertEqual(mcs.numAtoms, 3)
        self.assertEqual(mcs.numBonds, 3)

class TestCase(unittest.TestCase):

    def setUp(self):
        pass

    def test1(self):
        Common.test1(self)

    def test1PythonImpl(self):
        Common.test1(self,
            AtomTyper=CompareElements,
            BondTyper=CompareOrder)

    def test2(self):
        Common.test2(self)

    def test2PythonImpl(self):
        Common.test2(self,
            AtomTyper=CompareElements,
            BondTyper=CompareOrder)

    def test3IsotopeMatch(self):
        Common.test3IsotopeMatch(self)

    def test3IsotopeMatchPythonImpl(self):
        Common.test3IsotopeMatch(self,
            AtomTyper=CompareElements,
            BondTyper=CompareOrder)

    def test4RingMatches(self):
        Common.test4RingMatches(self)

    def test4RingMatchesPythonImpl(self):
        Common.test4RingMatches(self,
            AtomTyper=CompareElements,
            BondTyper=CompareOrder)

    def test5AnyMatch(self):
        Common.test5AnyMatch(self)

    def test5AnyMatchPythonImpl(self):
        Common.test5AnyMatch(self,
            AtomTyper=CompareElements,
            BondTyper=CompareOrder)

    def testAtomCompareAnyHeavyAtom(self):
        Common.testAtomCompareAnyHeavyAtom(self)

    def testAtomCompareAnyHeavyAtomPythonImpl(self):
        Common.testAtomCompareAnyHeavyAtom(self,
            AtomTyper=CompareElements,
            BondTyper=CompareOrder)

    def testAtomCompareAnyHeavyAtom1(self):
        Common.testAtomCompareAnyHeavyAtom1(self)

    def testAtomCompareAnyHeavyAtom1PythonImpl(self):
        Common.testAtomCompareAnyHeavyAtom1(self,
            AtomTyper=CompareElements,
            BondTyper=CompareOrder)

    def test6MatchValences(self):
        Common.test6MatchValences(self)

    def test6MatchValencesPythonImpl(self):
        Common.test6MatchValences(self,
            AtomTyper=CompareElements,
            BondTyper=CompareOrder)

    def test7Seed(self):
        Common.test7Seed(self)

    def test7SeedPythonImpl(self):
        Common.test7Seed(self,
            AtomTyper=CompareElements,
            BondTyper=CompareOrder)

    def test8MatchParams(self):
        Common.test8MatchParams(self)

    def test8MatchParamsPythonImpl(self):
        Common.test8MatchParams(self,
            AtomTyper=CompareElements,
            BondTyper=CompareOrder)

    def test9MatchCharge(self):
        Common.test9MatchCharge(self)

    def test9MatchChargePythonImpl(self):
        Common.test9MatchCharge(self,
            AtomTyper=CompareElements,
            BondTyper=CompareOrder)

    def test10MatchChargeAndParams(self):
        Common.test10MatchChargeAndParams(self)

    def test10MatchChargeAndParamsPythonImpl(self):
        Common.test10MatchChargeAndParams(self,
            AtomTyper=CompareElements,
            BondTyper=CompareOrder)

    def test11Github2034(self):
        Common.test11Github2034(self)

    def test11Github2034PythonImpl(self):
        Common.test11Github2034(self,
            AtomTyper=CompareElements,
            BondTyper=CompareOrder)

    def test12MCSAtomCompareExceptions(self):
        ps = rdFMCS.MCSParameters()
        smis = ['CCCCC', 'CCC1CCCCC1']
        ms = [Chem.MolFromSmiles(x) for x in smis]
        self.assertRaises(TypeError, lambda ps: setattr(ps, "AtomTyper",
                          AtomCompareCompareIsInt()))
        ps.AtomTyper = AtomCompareNoCompare()
        self.assertRaises(TypeError, lambda ms, ps: rdFMCS.FindMCS(ms, ps))

    def test13MCSAtomCompareUserData(self):
        smis = ['CCOCCOC', 'CCNCCCC']
        ms = [Chem.MolFromSmiles(x) for x in smis]
        ps = rdFMCS.MCSParameters()
        ps.AtomTyper = AtomCompareUserData()
        ps.AtomTyper.setMatchAnyHet(False)
        mcs = rdFMCS.FindMCS(ms, ps)
        self.assertEqual(mcs.numAtoms, 2)
        ps.AtomTyper.setMatchAnyHet(True)
        mcs = rdFMCS.FindMCS(ms, ps)
        self.assertEqual(mcs.numAtoms, 5)

    def test14MCSBondCompareExceptions(self):
        ps = rdFMCS.MCSParameters()
        smis = ['CCCCC', 'CCC1CCCCC1']
        ms = [Chem.MolFromSmiles(x) for x in smis]
        self.assertRaises(TypeError, lambda ps: setattr(ps, "BondTyper",
                          BondCompareCompareIsInt()))
        ps.BondTyper = BondCompareNoCompare()
        self.assertRaises(TypeError, lambda ms, ps: rdFMCS.FindMCS(ms, ps))

    def test15MCSBondCompareUserData(self):
        smis = ['C1CC=CCC1', 'c1ccccc1']
        ms = [Chem.MolFromSmiles(x) for x in smis]
        ps = rdFMCS.MCSParameters()
        ps.BondTyper = BondCompareUserData()
        ps.BondCompareParameters.CompleteRingsOnly = True
        ps.BondTyper.setIgnoreAromatization(False)
        mcs = rdFMCS.FindMCS(ms, ps)
        self.assertEqual(mcs.numAtoms, 0)
        ps.BondTyper.setIgnoreAromatization(True)
        mcs = rdFMCS.FindMCS(ms, ps)
        self.assertEqual(mcs.numAtoms, 6)

    def test16MCSProgressCallbackExceptions(self):
        ps = rdFMCS.MCSParameters()
        smis = ['CCCC(C)CC(CC)CC', 'OC(N)CC(C)CC(CC)CC']
        ms = [Chem.MolFromSmiles(x) for x in smis]
        self.assertRaises(TypeError, lambda ps: setattr(ps, "ProgressCallback",
                          ProgressCallbackCallbackIsInt()))
        ps.ProgressCallback = ProgressCallbackNoCallback()
        self.assertRaises(TypeError, lambda ms, ps: rdFMCS.FindMCS(ms, ps))

    def test17MCSProgressCallbackCancel(self):
        ps = rdFMCS.MCSParameters()
        smis = ['CCCC(C)CC(CC)CC', 'OC(N)CC(C)CC(CC)CC']
        ms = [Chem.MolFromSmiles(x) for x in smis]
        ps.AtomTyper = CompareElements()
        ps.ProgressCallback = ProgressCallback(self)
        mcs = rdFMCS.FindMCS(ms, ps)
        self.assertTrue(mcs.canceled)
        self.assertEqual(ps.ProgressCallback.callCount, 3)

if __name__ == "__main__":
    unittest.main()
