import unittest
import sys
from io import StringIO
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
    def __call__(self, p, mol1, atom1, mol2, atom2):
        if (p.MatchChiralTag and not self.CheckAtomChirality(p, mol1, atom1, mol2, atom2)):
            return False
        if (p.MatchFormalCharge and not self.CheckAtomCharge(p, mol1, atom1, mol2, atom2)):
            return False
        if (p.RingMatchesRingOnly):
            return self.CheckAtomRingMatch(p, mol1, atom1, mol2, atom2)
        return True

class CompareAnyHeavyAtom(CompareAny):
    def __call__(self, p, mol1, atom1, mol2, atom2):
        a1 = mol1.GetAtomWithIdx(atom1)
        a2 = mol2.GetAtomWithIdx(atom2)
        # Any atom, including H, matches another atom of the same type,  according to
        # the other flags
        if (a1.GetAtomicNum() == a2.GetAtomicNum() or
            (a1.GetAtomicNum() > 1 and a2.GetAtomicNum() > 1)):
            return CompareAny.__call__(self, p, mol1, atom1, mol2, atom2)
        return False

class CompareElements(rdFMCS.MCSAtomCompare):
    def __call__(self, p, mol1, atom1, mol2, atom2):
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
    def __call__(self, p, mol1, atom1, mol2, atom2):
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
    def __call__(self, p, mol1, bond1, mol2, bond2):
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
    __call__ = 1

class AtomCompareNoCompare(rdFMCS.MCSAtomCompare):
    pass

class AtomCompareUserData(rdFMCS.MCSAtomCompare):
    def __init__(self):
        super().__init__()
        self._matchAnyHet = False
    def setMatchAnyHet(self, v):
        self._matchAnyHet = v
    def __call__(self, p, mol1, atom1, mol2, atom2):
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
    __call__ = 1

class BondCompareNoCompare(rdFMCS.MCSBondCompare):
    pass

class BondCompareUserData(rdFMCS.MCSBondCompare):
    def __init__(self):
        super().__init__()
        self.match = None
    def setIgnoreAromatization(self, v):
        self.match = BondMatchOrderMatrix(v)
    def __call__(self, p, mol1, bond1, mol2, bond2):
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
    __call__ = 1

class ProgressCallbackNoCallback(rdFMCS.MCSProgress):
    pass

class ProgressCallback(rdFMCS.MCSProgress):
    def __init__(self, parent):
        super().__init__()
        self.parent = parent
        self.callCount = 0
    def __call__(self, stat, params):
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
        params = rdFMCS.MCSParameters()
        for kw in ("AtomTyper", "BondTyper"):
            v = kwargs.get(kw, None)
            if v is not None:
                v_instance = v()
                setattr(params, kw, v_instance)
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
        self.assertEqual(r.smartsString, "[#6]1-[#6]-[#6]-[#6]-1")
        self.assertEqual(r.numAtoms, 4)
        self.assertEqual(r.numBonds, 4)
        if kwargs:
            params = Common.getParams(**kwargs)
            params.InitialSeed = 'C1OC1'
            params.AtomCompareParameters.RingMatchesRingOnly = True
            params.BondCompareParameters.RingMatchesRingOnly = True
            r = rdFMCS.FindMCS(ms, params)
        else:
            r = rdFMCS.FindMCS(ms, seedSmarts='C1OC1', ringMatchesRingOnly=True)
        self.assertEqual(r.smartsString, "[#6&R]1-&@[#6&R]-&@[#6&R]-&@[#6&R]-&@1")
        self.assertEqual(r.numAtoms, 4)
        self.assertEqual(r.numBonds, 4)
        if kwargs:
            params = Common.getParams(**kwargs)
            params.InitialSeed = 'C1OC1'
            params.BondCompareParameters.CompleteRingsOnly = True
            r = rdFMCS.FindMCS(ms, params)
        else:
            r = rdFMCS.FindMCS(ms, seedSmarts='C1OC1', completeRingsOnly=True)
        self.assertEqual(r.smartsString, "[#6]1-&@[#6]-&@[#6]-&@[#6]-&@1")
        self.assertEqual(r.numAtoms, 4)
        self.assertEqual(r.numBonds, 4)

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
        ps.BondCompareParameters.CompleteRingsOnly = True
        ps.BondCompareParameters.MatchFusedRings = True
        mcs = rdFMCS.FindMCS(ms, ps)
        self.assertEqual(mcs.numAtoms, 3)

        ps = Common.getParams(**kwargs)
        ps.BondCompareParameters.CompleteRingsOnly = True
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

    def test19MCS3d(self, **kwargs):
        block1 = """
     RDKit          3D

 17 17  0  0  0  0  0  0  0  0999 V2000
    0.1592    0.8577    0.8639 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9090    0.9385   -0.2218 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1866   -0.4482   -0.7339 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0705   -1.1345   -1.1348 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8460   -1.2511   -0.0884 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3508    0.1733    0.1788 C   0  0  1  0  0  0  0  0  0  0  0  0
    1.6294    0.7269   -1.0491 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2088    0.1739    1.6399 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.4334    1.8588    1.2198 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8657    1.3295    0.1967 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5185    1.5822   -1.0340 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7155   -1.0290    0.0396 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8809   -0.3833   -1.6049 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.6828   -1.9159   -0.3875 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.3276   -1.6370    0.8016 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.2158    0.0912    0.8546 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.6002    1.6926   -1.1014 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  1  0
  4  5  1  0
  5  6  1  0
  6  7  1  0
  6  1  1  0
  1  8  1  0
  1  9  1  0
  2 10  1  0
  2 11  1  0
  3 12  1  0
  3 13  1  0
  5 14  1  0
  5 15  1  0
  6 16  1  1
  7 17  1  0
M  END


"""
        block2 = """
     RDKit          3D

 17 17  0  0  0  0  0  0  0  0999 V2000
    0.1592    0.8577    0.8639 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.9090    0.9385   -0.2218 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.1866   -0.4482   -0.7339 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.0705   -1.1345   -1.1348 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.8460   -1.2511   -0.0884 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.3508    0.1733    0.1788 C   0  0  2  0  0  0  0  0  0  0  0  0
    2.4771    0.1924    0.9509 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.1638    0.2755    1.7318 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.4849    1.8825    1.0902 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.4598    1.5686   -1.0290 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8337    1.3595    0.1902 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.8991   -0.3820   -1.5700 H   0  0  0  0  0  0  0  0  0  0  0  0
   -1.6394   -1.0265    0.1178 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.6816   -1.9167   -0.3664 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.3727   -1.5979    0.8445 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.4834    0.6899   -0.7827 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.2370    0.0336    1.8861 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  3  4  1  0
  4  5  1  0
  5  6  1  0
  6  7  1  0
  6  1  1  0
  1  8  1  0
  1  9  1  0
  2 10  1  0
  2 11  1  0
  3 12  1  0
  3 13  1  0
  5 14  1  0
  5 15  1  0
  6 16  1  6
  7 17  1  0
M  END


"""
        m1 = Chem.MolFromMolBlock(block1, removeHs=False)
        m2 = Chem.MolFromMolBlock(block2, removeHs=False)
        ps = Common.getParams(**kwargs)
        ps.AtomCompareParameters.MaxDistance = 1.0
        mcs = rdFMCS.FindMCS([m1, m2], ps)
        self.assertEqual(mcs.numAtoms, 14)
        self.assertEqual(mcs.numBonds, 14)


class TestCase(unittest.TestCase):

    def setUp(self):
        pass

    def test1(self):
        Common.test1(self)

    def test1PythonImpl(self):
        Common.test1(self,
            AtomTyper=CompareElements,
            BondTyper=CompareOrder)

    # DEPRECATED: remove from here in release 2021.01
    def test1PythonImplDeprecated(self):
        atom_call = CompareElements.__call__
        setattr(CompareElements, "compare", CompareElements.__call__)
        delattr(CompareElements, "__call__")
        bond_call = CompareOrder.__call__
        setattr(CompareOrder, "compare", CompareOrder.__call__)
        delattr(CompareOrder, "__call__")
        Common.test1(self,
            AtomTyper=CompareElements,
            BondTyper=CompareOrder)
        setattr(CompareElements, "__call__", atom_call)
        delattr(CompareElements, "compare")
        setattr(CompareOrder, "__call__", bond_call)
        delattr(CompareOrder, "compare")

    def test1PythonImplDeprecatedTypo(self):
        atom_call = CompareElements.__call__
        setattr(CompareElements, "comparx", CompareElements.__call__)
        delattr(CompareElements, "__call__")
        bond_call = CompareOrder.__call__
        setattr(CompareOrder, "comparx", CompareOrder.__call__)
        delattr(CompareOrder, "__call__")
        self.assertRaises(TypeError, lambda self: Common.test1(self,
            AtomTyper=CompareElements,
            BondTyper=CompareOrder))
        setattr(CompareElements, "__call__", atom_call)
        delattr(CompareElements, "comparx")
        setattr(CompareOrder, "__call__", bond_call)
        delattr(CompareOrder, "comparx")
    # DEPRECATED: remove until here in release 2021.01

    def test2(self):
        Common.test2(self)

    def test2PythonImpl(self):
        Common.test2(self,
            AtomTyper=CompareElements,
            BondTyper=CompareOrder)

    def test2PythonImplAtomTyperOnly(self):
        Common.test2(self,
            AtomTyper=CompareElements)

    def test2PythonImplBondTyperOnly(self):
        Common.test2(self,
            BondTyper=CompareOrder)

    def test3IsotopeMatch(self):
        Common.test3IsotopeMatch(self)

    def test3IsotopeMatchPythonImpl(self):
        Common.test3IsotopeMatch(self,
            AtomTyper=CompareElements,
            BondTyper=CompareOrder)

    def test3IsotopeMatchPythonImplAtomTyperOnly(self):
        Common.test3IsotopeMatch(self,
            AtomTyper=CompareElements)

    def test3IsotopeMatchPythonImplBondTyperOnly(self):
        Common.test3IsotopeMatch(self,
            BondTyper=CompareOrder)

    def test4RingMatches(self):
        Common.test4RingMatches(self)

    def test4RingMatchesPythonImpl(self):
        Common.test4RingMatches(self,
            AtomTyper=CompareElements,
            BondTyper=CompareOrder)

    def test4RingMatchesPythonImplAtomTyperOnly(self):
        Common.test4RingMatches(self,
            AtomTyper=CompareElements)

    def test4RingMatchesPythonImplBondTyperOnly(self):
        Common.test4RingMatches(self,
            BondTyper=CompareOrder)

    def test5AnyMatch(self):
        Common.test5AnyMatch(self)

    def test5AnyMatchPythonImpl(self):
        Common.test5AnyMatch(self,
            AtomTyper=CompareElements,
            BondTyper=CompareOrder)

    def test5AnyMatchPythonImplAtomTyperOnly(self):
        Common.test5AnyMatch(self,
            AtomTyper=CompareElements)

    def test5AnyMatchPythonImplBondTyperOnly(self):
        Common.test5AnyMatch(self,
            BondTyper=CompareOrder)

    def testAtomCompareAnyHeavyAtom(self):
        Common.testAtomCompareAnyHeavyAtom(self)

    def testAtomCompareAnyHeavyAtomPythonImpl(self):
        Common.testAtomCompareAnyHeavyAtom(self,
            AtomTyper=CompareElements,
            BondTyper=CompareOrder)

    def testAtomCompareAnyHeavyAtomPythonImplAtomTyperOnly(self):
        Common.testAtomCompareAnyHeavyAtom(self,
            AtomTyper=CompareElements)

    def testAtomCompareAnyHeavyAtomPythonImplBondTyperOnly(self):
        Common.testAtomCompareAnyHeavyAtom(self,
            BondTyper=CompareOrder)

    def testAtomCompareAnyHeavyAtom1(self):
        Common.testAtomCompareAnyHeavyAtom1(self)

    def testAtomCompareAnyHeavyAtom1PythonImpl(self):
        Common.testAtomCompareAnyHeavyAtom1(self,
            AtomTyper=CompareElements,
            BondTyper=CompareOrder)

    def testAtomCompareAnyHeavyAtom1PythonImplAtomTyperOnly(self):
        Common.testAtomCompareAnyHeavyAtom1(self,
            AtomTyper=CompareElements)

    def testAtomCompareAnyHeavyAtom1PythonImplBondTyperOnly(self):
        Common.testAtomCompareAnyHeavyAtom1(self,
            BondTyper=CompareOrder)

    def test6MatchValences(self):
        Common.test6MatchValences(self)

    def test6MatchValencesPythonImpl(self):
        Common.test6MatchValences(self,
            AtomTyper=CompareElements,
            BondTyper=CompareOrder)

    def test6MatchValencesPythonImplAtomTyperOnly(self):
        Common.test6MatchValences(self,
            AtomTyper=CompareElements)

    def test6MatchValencesPythonImplBondTyperOnly(self):
        Common.test6MatchValences(self,
            BondTyper=CompareOrder)

    def test7Seed(self):
        Common.test7Seed(self)

    def test7SeedPythonImpl(self):
        Common.test7Seed(self,
            AtomTyper=CompareElements,
            BondTyper=CompareOrder)

    def test7SeedPythonImplAtomTyperOnly(self):
        Common.test7Seed(self,
            AtomTyper=CompareElements)

    def test7SeedPythonImplBondTyperOnly(self):
        Common.test7Seed(self,
            BondTyper=CompareOrder)

    def test8MatchParams(self):
        Common.test8MatchParams(self)

    def test8MatchParamsPythonImpl(self):
        Common.test8MatchParams(self,
            AtomTyper=CompareElements,
            BondTyper=CompareOrder)

    def test8MatchParamsPythonImplAtomTyperOnly(self):
        Common.test8MatchParams(self,
            AtomTyper=CompareElements)

    def test8MatchParamsPythonImplBondTyperOnly(self):
        Common.test8MatchParams(self,
            BondTyper=CompareOrder)

    def test9MatchCharge(self):
        Common.test9MatchCharge(self)

    def test9MatchChargePythonImpl(self):
        Common.test9MatchCharge(self,
            AtomTyper=CompareElements,
            BondTyper=CompareOrder)

    def test9MatchChargePythonImplAtomTyperOnly(self):
        Common.test9MatchCharge(self,
            AtomTyper=CompareElements)

    def test9MatchChargePythonImplBondTyperOnly(self):
        Common.test9MatchCharge(self,
            BondTyper=CompareOrder)

    def test10MatchChargeAndParams(self):
        Common.test10MatchChargeAndParams(self)

    def test10MatchChargeAndParamsPythonImpl(self):
        Common.test10MatchChargeAndParams(self,
            AtomTyper=CompareElements,
            BondTyper=CompareOrder)

    def test10MatchChargeAndParamsPythonImplAtomTyperOnly(self):
        Common.test10MatchChargeAndParams(self,
            AtomTyper=CompareElements)

    def test10MatchChargeAndParamsPythonImplBondTyperOnly(self):
        Common.test10MatchChargeAndParams(self,
            BondTyper=CompareOrder)

    def test11Github2034(self):
        Common.test11Github2034(self)

    def test11Github2034PythonImpl(self):
        Common.test11Github2034(self,
            AtomTyper=CompareElements,
            BondTyper=CompareOrder)

    def test11Github2034PythonImplAtomTyperOnly(self):
        Common.test11Github2034(self,
            AtomTyper=CompareElements)

    def test11Github2034PythonImplBondTyperOnly(self):
        Common.test11Github2034(self,
            BondTyper=CompareOrder)

    def test12MCSAtomCompareExceptions(self):
        ps = rdFMCS.MCSParameters()
        smis = ['CCCCC', 'CCC1CCCCC1']
        ms = [Chem.MolFromSmiles(x) for x in smis]
        self.assertRaises(TypeError, lambda ps: setattr(ps, "AtomTyper",
                          AtomCompareCompareIsInt()))
        self.assertRaises(TypeError, lambda ps: setattr(ps, "AtomTyper",
                          AtomCompareNoCompare()))

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
        self.assertRaises(TypeError, lambda ps: setattr(ps, "BondTyper",
                          BondCompareNoCompare()))

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
        self.assertRaises(TypeError, lambda ps: setattr(ps, "ProgressCallback",
                          ProgressCallbackNoCallback()))

    def test17MCSProgressCallbackCancel(self):
        ps = rdFMCS.MCSParameters()
        smis = ['CCCC(C)CC(CC)CC', 'OC(N)CC(C)CC(CC)CC']
        ms = [Chem.MolFromSmiles(x) for x in smis]
        ps.AtomTyper = CompareElements()
        ps.ProgressCallback = ProgressCallback(self)
        mcs = rdFMCS.FindMCS(ms, ps)
        self.assertTrue(mcs.canceled)
        self.assertEqual(ps.ProgressCallback.callCount, 3)

    # DEPRECATED: remove from here in release 2021.01
    def test17MCSProgressCallbackCancelDeprecated(self):
        callback = ProgressCallback.__call__
        setattr(ProgressCallback, "callback", ProgressCallback.__call__)
        delattr(ProgressCallback, "__call__")
        self.test17MCSProgressCallbackCancel()
        setattr(ProgressCallback, "__call__", callback)
        delattr(ProgressCallback, "callback")

    def test17MCSProgressCallbackCancelDeprecatedTypo(self):
        callback = ProgressCallback.__call__
        setattr(ProgressCallback, "callbacx", ProgressCallback.__call__)
        delattr(ProgressCallback, "__call__")
        self.assertRaises(TypeError, self.test17MCSProgressCallbackCancel)
        setattr(ProgressCallback, "__call__", callback)
        delattr(ProgressCallback, "callbacx")
    # DEPRECATED: remove until here in release 2021.01

    def test18GitHub3693(self):
        mols = [Chem.MolFromSmiles(smi) for smi in [
                "Nc1ccc(O)cc1c1ccc2ccccc2c1", "Oc1cnc(NC2CCC2)c(c1)c1ccc2ccccc2c1"]]
        params = rdFMCS.MCSParameters()
        res = rdFMCS.FindMCS(mols, params)
        self.assertEqual(res.numAtoms, 17)
        self.assertEqual(res.numBonds, 18)
        self.assertEqual(res.smartsString, "[#7]-,:[#6]:[#6](:[#6]:[#6](:[#6])-[#8])-[#6]1:[#6]:[#6]:[#6]2:[#6](:[#6]:1):[#6]:[#6]:[#6]:[#6]:2")

        params = rdFMCS.MCSParameters()
        params.BondCompareParameters.CompleteRingsOnly = True
        res = rdFMCS.FindMCS(mols, params)
        self.assertEqual(res.numAtoms, 11)
        self.assertEqual(res.numBonds, 12)
        self.assertEqual(res.smartsString, "[#6]-&!@[#6]1:&@[#6]:&@[#6]:&@[#6]2:&@[#6](:&@[#6]:&@1):&@[#6]:&@[#6]:&@[#6]:&@[#6]:&@2")

        params = rdFMCS.MCSParameters()
        params.AtomCompareParameters.CompleteRingsOnly = True
        params.BondCompareParameters.CompleteRingsOnly = True
        res = rdFMCS.FindMCS(mols, params)
        self.assertEqual(res.numAtoms, 10)
        self.assertEqual(res.numBonds, 11)
        self.assertEqual(res.smartsString, "[#6&R]1:&@[#6&R]:&@[#6&R]:&@[#6&R]2:&@[#6&R](:&@[#6&R]:&@1):&@[#6&R]:&@[#6&R]:&@[#6&R]:&@[#6&R]:&@2")

        params = rdFMCS.MCSParameters()
        params.AtomCompareParameters.CompleteRingsOnly = True
        # this will automatically be set to True
        params.BondCompareParameters.CompleteRingsOnly = False
        res = rdFMCS.FindMCS(mols, params)
        self.assertEqual(res.numAtoms, 10)
        self.assertEqual(res.numBonds, 11)
        self.assertEqual(res.smartsString, "[#6&R]1:&@[#6&R]:&@[#6&R]:&@[#6&R]2:&@[#6&R](:&@[#6&R]:&@1):&@[#6&R]:&@[#6&R]:&@[#6&R]:&@[#6&R]:&@2")

    def test19MCS3d(self):
        Common.test19MCS3d(self)

    def test20AtomCompareCompleteRingsOnly(self):
        mols = [Chem.MolFromSmiles(smi) for smi in ["C1CCCC1C", "C1CCCC1C1CCCCC1"]]
        params = rdFMCS.MCSParameters()
        params.AtomCompareParameters.CompleteRingsOnly = True
        res = rdFMCS.FindMCS(mols, params)
        self.assertEqual(res.numAtoms, 5)
        self.assertEqual(res.numBonds, 5)
        self.assertEqual(res.smartsString, "[#6&R]1-&@[#6&R]-&@[#6&R]-&@[#6&R]-&@[#6&R]-&@1")

        params = rdFMCS.MCSParameters()
        params.AtomCompareParameters.CompleteRingsOnly = True
        # this will automatically be set to True
        params.BondCompareParameters.CompleteRingsOnly = False
        res = rdFMCS.FindMCS(mols, params)
        self.assertEqual(res.numAtoms, 5)
        self.assertEqual(res.numBonds, 5)
        self.assertEqual(res.smartsString, "[#6&R]1-&@[#6&R]-&@[#6&R]-&@[#6&R]-&@[#6&R]-&@1")

if __name__ == "__main__":
    unittest.main()
