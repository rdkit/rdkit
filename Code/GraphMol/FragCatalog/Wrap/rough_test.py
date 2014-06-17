# $Id$
#
#  Copyright (C) 2003-2005  Rational Discovery LLC
#         All Rights Reserved
#
""" This is a rough coverage test of the python wrapper

it's intended to be shallow, but broad

"""
import unittest,os
from rdkit.six.moves import cPickle
from rdkit import RDConfig
from rdkit.RDLogger import logger
logger=logger()
from rdkit import Chem
from rdkit.Chem import FragmentCatalog
from rdkit import DataStructs


class TestCase(unittest.TestCase):

    def setUp(self) :
        self.fName=os.path.join(RDConfig.RDBaseDir,'Code','GraphMol',
                                'FragCatalog','test_data','funcGroups.txt')
        self.smiName=os.path.join(RDConfig.RDBaseDir,'Code','GraphMol',
                                  'FragCatalog','test_data','mols.smi')

    def test0Params(self) :
        fparams = FragmentCatalog.FragCatParams(1, 6, self.fName, 1.0e-8)

        ctype = fparams.GetTypeString()
        assert(ctype == "Fragment Catalog Parameters") 
        assert(fparams.GetLowerFragLength() == 1)
        assert(fparams.GetUpperFragLength() == 6)
        
        ngps = fparams.GetNumFuncGroups()
        assert ngps==15
        for i in range(ngps) :
          mol = fparams.GetFuncGroup(i)
            
    def test1Catalog(self) :
        fparams = FragmentCatalog.FragCatParams(1, 6, self.fName, 1.0e-8)
        fcat = FragmentCatalog.FragCatalog(fparams)
        assert(fcat.GetNumEntries() == 0)
        assert( fcat.GetFPLength() == 0)
        
        nparams = fcat.GetCatalogParams()
        assert(nparams.GetLowerFragLength() == 1)
        assert(nparams.GetUpperFragLength() == 6)
        
    def test2Generator(self) :
        fparams = FragmentCatalog.FragCatParams(1, 6, self.fName, 1.0e-8)
        fcat = FragmentCatalog.FragCatalog(fparams)
        fgen = FragmentCatalog.FragCatGenerator()
        suppl = Chem.SmilesMolSupplier(self.smiName," ",0,1,0)
        for mol in suppl:
            nent = fgen.AddFragsFromMol(mol, fcat)
        assert fcat.GetNumEntries()==21
        assert fcat.GetFPLength()==21
        for id in range(fcat.GetNumEntries()):
            assert fcat.GetEntryBitId(id)==id
            assert fcat.GetEntryOrder(id)==fcat.GetBitOrder(id)
            assert fcat.GetEntryDescription(id)==fcat.GetBitDescription(id)
            assert tuple(fcat.GetEntryFuncGroupIds(id))==tuple(fcat.GetBitFuncGroupIds(id))
            
    def test3FPgenerator(self) :
        with open(self.smiName,'r') as smiF:
            smiLines = smiF.readlines()
        fparams = FragmentCatalog.FragCatParams(1, 6, self.fName)
        fcat = FragmentCatalog.FragCatalog(fparams)
        fgen = FragmentCatalog.FragCatGenerator()
        suppl = Chem.SmilesMolSupplier(self.smiName," ",0,1,0)
        smiles = []
        for mol in suppl:
            nent = fgen.AddFragsFromMol(mol, fcat)
            smiles.append(Chem.MolToSmiles(mol))
        assert fcat.GetNumEntries()==21
        assert fcat.GetFPLength()==21,fcat.GetFPLength()        
        fpgen = FragmentCatalog.FragFPGenerator()
        obits = [3,2,3,3,2,3,5,5,5,4,5,6]
        obls = [(0,1,2),(1,3),(1,4,5),(1,6,7),(0,8),(0,6,9),(0,1,2,3,10),
                (0,1,2,8,11),(1,3,4,5,12),(1,4,5,13),(1,3,6,7,14),(0,1,6,7,9,15)]
        for i in range(len(smiles)):
            smi = smiles[i]
            mol = Chem.MolFromSmiles(smi)
            fp = fpgen.GetFPForMol(mol, fcat)
            if i < len(obits):
                assert fp.GetNumOnBits()==obits[i],'%s: %s'%(smi,str(fp.GetOnBits()))
            obl = fp.GetOnBits()
            if i < len(obls):
                assert tuple(obl)==obls[i],'%s: %s'%(smi,obl)

                
    def test4Serialize(self) :
        with open(self.smiName,'r') as smiF:
            smiLines = smiF.readlines()
        fparams = FragmentCatalog.FragCatParams(1, 6, self.fName)
        fcat = FragmentCatalog.FragCatalog(fparams)
        fgen = FragmentCatalog.FragCatGenerator()
        suppl = Chem.SmilesMolSupplier(self.smiName," ",0,1,0)
        smiles = []
        for mol in suppl:
            nent = fgen.AddFragsFromMol(mol, fcat)
            smiles.append(Chem.MolToSmiles(mol))
        assert fcat.GetNumEntries()==21
        assert fcat.GetFPLength()==21,fcat.GetFPLength()        
        pkl = cPickle.dumps(fcat)
        fcat2 = cPickle.loads(pkl)
        assert fcat2.GetNumEntries()==21
        assert fcat2.GetFPLength()==21,fcat2.GetFPLength()        
        fpgen = FragmentCatalog.FragFPGenerator()
        for i in range(len(smiles)):
            smi = smiles[i]
            mol = Chem.MolFromSmiles(smi)
            fp1 = fpgen.GetFPForMol(mol, fcat)
            fp2 = fpgen.GetFPForMol(mol, fcat2)
            assert fp1.GetNumOnBits()==fp2.GetNumOnBits()
            obl1 = fp1.GetOnBits()
            obl2 = fp2.GetOnBits()
            assert tuple(obl1)==tuple(obl2)

                
    def test5FPsize(self) :
        with open(self.smiName,'r') as smiF:
            smiLines = smiF.readlines()
        fparams = FragmentCatalog.FragCatParams(6, 6, self.fName)
        fcat = FragmentCatalog.FragCatalog(fparams)
        fgen = FragmentCatalog.FragCatGenerator()
        suppl = [Chem.MolFromSmiles('C1CCCOC1O')]
        for mol in suppl:
            nent = fgen.AddFragsFromMol(mol, fcat)
        assert fcat.GetFPLength()==1
        for i in range(fcat.GetFPLength()):
            assert fcat.GetBitOrder(i)==6
            assert fcat.GetBitDescription(i)=="C1CCOC<-O>C1",fcat.GetBitDescription(i)
            assert tuple(fcat.GetBitFuncGroupIds(i))==(1,)
            
        
    def test6DownEntries(self) :
        fparams = FragmentCatalog.FragCatParams(1, 6, self.fName, 1.0e-8)
        fcat = FragmentCatalog.FragCatalog(fparams)
        fgen = FragmentCatalog.FragCatGenerator()
        suppl = Chem.SmilesMolSupplier(self.smiName," ",0,1,0)
        for mol in suppl:
            nent = fgen.AddFragsFromMol(mol, fcat)
        assert fcat.GetNumEntries()==21
        assert fcat.GetFPLength()==21
        assert tuple(fcat.GetEntryDownIds(0))==(2,8,9,16)
        assert tuple(fcat.GetEntryDownIds(1))==(2,3,5,7)

    def test7Issue116(self):
        smiList = ['Cc1ccccc1']
        suppl = Chem.SmilesMolSupplierFromText('\n'.join(smiList),
                                               ',',0,-1,0)
        fparams = FragmentCatalog.FragCatParams(2, 2, self.fName, 1.0e-8)
        cat = FragmentCatalog.FragCatalog(fparams)
        fgen = FragmentCatalog.FragCatGenerator()
        for mol in suppl:
            nent = fgen.AddFragsFromMol(mol, cat)
        assert cat.GetFPLength()==2
        assert cat.GetBitDescription(0)=='ccC'
        fpgen = FragmentCatalog.FragFPGenerator()
        mol = Chem.MolFromSmiles('Cc1ccccc1')
        fp = fpgen.GetFPForMol(mol,cat)
        assert fp[0]
        assert fp[1]
        
        mol = Chem.MolFromSmiles('c1ccccc1-c1ccccc1')
        fp = fpgen.GetFPForMol(mol,cat)
        assert not fp[0]
        assert fp[1]

    def test8Issue118(self):
        smiList = ['CCN(C(N)=O)N=O']
        fName=os.path.join(RDConfig.RDDataDir,'FunctionalGroups.txt')
        suppl = Chem.SmilesMolSupplierFromText('\n'.join(smiList),
                                               ',',0,-1,0)
        fparams = FragmentCatalog.FragCatParams(2, 4, fName, 1.0e-8)
        cat = FragmentCatalog.FragCatalog(fparams)
        fgen = FragmentCatalog.FragCatGenerator()
        for mol in suppl:
            nent = fgen.AddFragsFromMol(mol, cat)
        assert cat.GetFPLength()==1
        assert cat.GetBitDescription(0)=='CCN(<-C(=O)N>)<-N=O>'


if __name__ == '__main__':
    unittest.main()
