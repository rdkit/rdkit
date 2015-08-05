from rdkit import RDConfig
import unittest
from rdkit import Chem
from rdkit.Chem import rdMMPA

def natoms(tpl):
    return tuple(x.GetNumAtoms() if x is not None else 0 for x in tpl)

class TestCase(unittest.TestCase):
    def setUp(self) :
        pass

    def test1(self) :
        m = Chem.MolFromSmiles('c1ccccc1OC')
        frags = rdMMPA.FragmentMol(m)
        self.failUnlessEqual(len(frags),3)
        for frag in frags:
            self.failUnlessEqual(len(frag),2)        
        frags = sorted(frags,key=natoms)
        self.failUnlessEqual(frags[0][0],None)
        self.failUnlessEqual(frags[1][0],None)
        self.failIfEqual(frags[2][0],None)
        self.failIfEqual(frags[0][1],None)
        self.failIfEqual(frags[1][1],None)
        self.failIfEqual(frags[2][1],None)

        self.failUnlessEqual(frags[0][1].GetNumAtoms(),m.GetNumAtoms()+2)
        self.failUnlessEqual(frags[1][1].GetNumAtoms(),m.GetNumAtoms()+2)

        fs = Chem.GetMolFrags(frags[0][1],asMols=True)
        self.failUnlessEqual(len(fs),2)
        self.failUnlessEqual(Chem.MolToSmiles(fs[0],True),'c1ccc([*:1])cc1')
        self.failUnlessEqual(Chem.MolToSmiles(fs[1],True),'CO[*:1]')

        fs = Chem.GetMolFrags(frags[1][1],asMols=True)
        self.failUnlessEqual(len(fs),2)
        self.failUnlessEqual(Chem.MolToSmiles(fs[0],True),'c1ccc(O[*:1])cc1')
        self.failUnlessEqual(Chem.MolToSmiles(fs[1],True),'C[*:1]')

        fs = Chem.GetMolFrags(frags[2][0],asMols=True)
        self.failUnlessEqual(len(fs),1)
        self.failUnlessEqual(Chem.MolToSmiles(fs[0],True),'O([*:1])[*:2]')
        fs = Chem.GetMolFrags(frags[2][1],asMols=True)
        self.failUnlessEqual(len(fs),2)
        self.failUnlessEqual(Chem.MolToSmiles(fs[0],True),'c1ccc([*:1])cc1')
        self.failUnlessEqual(Chem.MolToSmiles(fs[1],True),'C[*:2]')

        
    def test2(self) :
        m = Chem.MolFromSmiles('c1ccccc1OC')
        frags = rdMMPA.FragmentMol(m,resultsAsMols=False)
        self.failUnlessEqual(len(frags),3)
        for frag in frags:
            self.failUnlessEqual(len(frag),2)        
        frags = sorted(frags)
        self.failUnlessEqual(frags[0][0],'')
        self.failUnlessEqual(frags[1][0],'')
        self.failIfEqual(frags[2][0],'')
        self.failIfEqual(frags[0][1],'')
        self.failIfEqual(frags[1][1],'')
        self.failIfEqual(frags[2][1],'')
        
        self.failUnlessEqual(frags[0][1],'CO[*:1].c1ccc(cc1)[*:1]')
        self.failUnlessEqual(frags[1][1],'C[*:1].c1ccc(cc1)O[*:1]')
        self.failUnlessEqual(frags[2][0],'O([*:1])[*:2]')
        self.failUnlessEqual(frags[2][1],'C[*:2].c1ccc([*:1])cc1')

    def test3(self) :
        m = Chem.MolFromSmiles('c1ccccc1OC')
        frags = rdMMPA.FragmentMol(m,resultsAsMols=False,pattern='cO')
        self.failUnlessEqual(len(frags),1)
        for frag in frags:
            self.failUnlessEqual(len(frag),2)        
        frags = sorted(frags)
        self.failUnlessEqual(frags[0][0],'')
        self.failIfEqual(frags[0][1],'')
        
        self.failUnlessEqual(frags[0][1],'CO[*:1].c1ccc(cc1)[*:1]')

    def test4(self): # currently failing
        m = Chem.MolFromSmiles('Cc1ccccc1NC(=O)C(C)[NH+]1CCCC1') # ZINC00000051
        frags = rdMMPA.FragmentMol(m,resultsAsMols=False)
        for frag in sorted(frags):
            print(frag)
        cores = set(x[0] for x in frags)
        self.failUnless('C([*:1])([*:2])[*:3]' in cores)
        self.failUnless('O=C(N[*:3])C([*:1])[*:2]' in cores)
        self.failUnlessEqual(len(frags),18)
        for frag in frags:
            self.failUnlessEqual(len(frag),2)        



if __name__=="__main__":
    unittest.main()
