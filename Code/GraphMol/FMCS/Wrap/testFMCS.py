from rdkit import RDConfig
import os,sys
import unittest
from rdkit import Chem
from rdkit.Chem import rdFMCS


class TestCase(unittest.TestCase):
    def setUp(self) :
        pass
    def test1(self):
        smis=(    "Cc1nc(CN(C(C)c2ncccc2)CCCCN)ccc1 CHEMBL1682991", #-- QUERY
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
        mcs = rdFMCS.FindMCS(ms)
        self.assertEqual(mcs.numBonds,21)
        self.assertEqual(mcs.numAtoms,21)
        self.assertEqual(mcs.smartsString,'[#6](:[#6]:[#6]):[#6]:[#7]:[#6]-[#6]-[#7](-[#6](-[#6])-[#6]1:[#6]:[#6]:[#6]:[#6]:[#7]:1)-[#6]-[#6]-[#6]-[#6]-[#7]')
        qm = Chem.MolFromSmarts(mcs.smartsString)
        self.failUnless(qm is not None)
        for m in ms:
            self.failUnless(m.HasSubstructMatch(qm))
        
        mcs = rdFMCS.FindMCS(ms,threshold=0.8)
        self.assertEqual(mcs.numBonds,21)
        self.assertEqual(mcs.numAtoms,21)
        self.assertEqual(mcs.smartsString,'[#6](:[#6]:[#6]):[#6]:[#7]:[#6]-[#6]-[#7](-[#6](-[#6])-[#6]1:[#6]:[#6]:[#6]:[#6]:[#7]:1)-[#6]-[#6]-[#6]-[#6]-[#7]')
        qm = Chem.MolFromSmarts(mcs.smartsString)
        self.failUnless(qm is not None)
        for m in ms:
            self.failUnless(m.HasSubstructMatch(qm))

    def test2(self):
        smis=(
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
        mcs = rdFMCS.FindMCS(ms)
        self.assertEqual(mcs.numBonds,9)
        self.assertEqual(mcs.numAtoms,10)
        qm = Chem.MolFromSmarts(mcs.smartsString)
        self.failUnless(qm is not None)
        for m in ms:
            self.failUnless(m.HasSubstructMatch(qm))
        # smarts too hard to canonicalize this
        #self.assertEqual(mcs.smartsString,'[#6]-,:[#6]-,:[#6]-,:[#6]-,:[#6](-[#6]-[#8]-[#6]:,-[#6])-,:[#6]')

        mcs = rdFMCS.FindMCS(ms,threshold=0.8)
        self.assertEqual(mcs.numBonds,20)
        self.assertEqual(mcs.numAtoms,19)
        qm = Chem.MolFromSmarts(mcs.smartsString)
        self.failUnless(qm is not None)
        nHits=0
        for m in ms:
            if m.HasSubstructMatch(qm):
                nHits+=1
        self.failUnless(nHits>=int(0.8*len(smis)))
        # smarts too hard to canonicalize this
        #self.assertEqual(mcs.smartsString,'[#6]1:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#6](-[#8]-[#6]-[#6]-[#7]-[#6]-[#6])-[#6]2:[#6]:[#6]:[#6]:[#6]:[#6]:2')

    def test3IsotopeMatch(self):
        smis=(
            "CC[14NH2]",
            "CC[14CH3]",
            )
        
        ms = [Chem.MolFromSmiles(x) for x in smis]
        mcs = rdFMCS.FindMCS(ms)
        self.assertEqual(mcs.numBonds,1)
        self.assertEqual(mcs.numAtoms,2)
        qm = Chem.MolFromSmarts(mcs.smartsString)

        mcs = rdFMCS.FindMCS(ms,atomCompare=rdFMCS.AtomCompare.CompareIsotopes)
        self.assertEqual(mcs.numBonds,2)
        self.assertEqual(mcs.numAtoms,3)
        qm = Chem.MolFromSmarts(mcs.smartsString)

        self.failUnless(Chem.MolFromSmiles('CC[14CH3]').HasSubstructMatch(qm))
        self.failIf(Chem.MolFromSmiles('CC[13CH3]').HasSubstructMatch(qm))
        self.failUnless(Chem.MolFromSmiles('OO[14CH3]').HasSubstructMatch(qm))
        self.failIf(Chem.MolFromSmiles('O[13CH2][14CH3]').HasSubstructMatch(qm))

    def test4RingMatches(self):
        smis = ['CCCCC','CCC1CCCCC1']
        ms = [Chem.MolFromSmiles(x) for x in smis]
        mcs = rdFMCS.FindMCS(ms)
        self.assertEqual(mcs.numBonds,4)
        self.assertEqual(mcs.numAtoms,5)
        self.assertEqual(mcs.smartsString,'[#6]-[#6]-[#6]-[#6]-[#6]')

        mcs = rdFMCS.FindMCS(ms,completeRingsOnly=True)
        self.assertEqual(mcs.numBonds,2)
        self.assertEqual(mcs.numAtoms,3)
        self.assertEqual(mcs.smartsString,'[#6]-[#6]-[#6]')

        mcs = rdFMCS.FindMCS(ms,ringMatchesRingOnly=True)
        self.assertEqual(mcs.numBonds,2)
        self.assertEqual(mcs.numAtoms,3)
        self.assertEqual(mcs.smartsString,'[#6]-[#6]-[#6]')

        smis = ['CC1CCC1','CCC1CCCCC1']
        ms = [Chem.MolFromSmiles(x) for x in smis]
        mcs = rdFMCS.FindMCS(ms)
        self.assertEqual(mcs.numBonds,4)
        self.assertEqual(mcs.numAtoms,5)
        self.assertEqual(mcs.smartsString,'[#6]-[#6](-[#6]-[#6])-[#6]')

        mcs = rdFMCS.FindMCS(ms,completeRingsOnly=True)
        self.assertEqual(mcs.numBonds,1)
        self.assertEqual(mcs.numAtoms,2)
        self.assertEqual(mcs.smartsString,'[#6]-[#6]')

        mcs = rdFMCS.FindMCS(ms,ringMatchesRingOnly=True)
        self.assertEqual(mcs.numBonds,4)
        self.assertEqual(mcs.numAtoms,5)
        self.assertEqual(mcs.smartsString,'[#6]-[#6](-[#6]-[#6])-[#6]')

    def test5AnyMatch(self):
        smis=('c1ccccc1C',
              'c1ccccc1O',
              'c1ccccc1Cl'
            )
        ms = [Chem.MolFromSmiles(x) for x in smis]
        mcs = rdFMCS.FindMCS(ms,atomCompare=rdFMCS.AtomCompare.CompareAny)
        self.assertEqual(mcs.numBonds,7)
        self.assertEqual(mcs.numAtoms,7)
        qm = Chem.MolFromSmarts(mcs.smartsString)

        for m in ms:
            self.assertTrue(m.HasSubstructMatch(qm))


        smis=('c1cccnc1C',
              'c1cnncc1O',
              'c1cccnc1Cl'
            )
        ms = [Chem.MolFromSmiles(x) for x in smis]
        mcs = rdFMCS.FindMCS(ms,atomCompare=rdFMCS.AtomCompare.CompareAny)
        self.assertEqual(mcs.numBonds,7)
        self.assertEqual(mcs.numAtoms,7)
        qm = Chem.MolFromSmarts(mcs.smartsString)

        for m in ms:
            self.assertTrue(m.HasSubstructMatch(qm))

    def test6MatchValences(self):
        ms = (Chem.MolFromSmiles('NC1OC1'),Chem.MolFromSmiles('C1OC1[N+](=O)[O-]'))
        mcs = rdFMCS.FindMCS(ms)
        self.assertEqual(mcs.numBonds,4)
        self.assertEqual(mcs.numAtoms,4)
        mcs = rdFMCS.FindMCS(ms,matchValences=True)
        self.assertEqual(mcs.numBonds,3)
        self.assertEqual(mcs.numAtoms,3)


        
        
            
if __name__=="__main__":
    unittest.main()
