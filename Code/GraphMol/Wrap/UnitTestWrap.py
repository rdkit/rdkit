# Copyright (C) 2001 greg Landrum and Rational Discovery LLC

"""basic unit testing code for the RDKit Wrapper

"""
import RDConfig
import unittest,os,sys,re
from Chem import rdchem as rdmol
from Chem.rdchem import Atom,Bond,Mol
from Numeric import *
import cPickle

def Smi2CDXML(smi):
  mol = rdmol.MolFromSmiles(smi)
  d = rdmol.MolToCDXML(mol)
  return d

_FTOL=1e-4
def _fmatcmp(m1,m2):
  return max(max(abs(m1-m2)))

class TestCase(unittest.TestCase):
  def setUp(self):
    #print '\n%s: '%self.shortDescription()
    self.utilsPath = '%s/Code/GraphMol/Utils'%(RDConfig.RDBaseDir)
    self.substructPath = '%s/Code/GraphMol/Substruct'%(RDConfig.RDBaseDir)
    self.dataPath = '%s/Code/GraphMol/Wrap/test_data'%(RDConfig.RDBaseDir)
  def _bulkSmilesTest(self,lines):
    splitExpr = re.compile('[\t\ ]')
    for line in lines:
      smi = splitExpr.split(line)[1]
      if len(smi):
        if smi[-1] == '\n': smi = smi[:-1]
        if smi[-1] == '\r': smi = smi[:-1]
        ok = 1
        try:
          m = rdmol.MolFromSmiles(smi)
        except:
          ok = 0
        assert ok,"parse of smiles %s failed"%(smi)

        canonSmi = rdmol.MolToSmiles(m)
        m = rdmol.MolFromSmiles(canonSmi)
        try:
          newCanon = rdmol.MolToSmiles(m)
        except:
          ok = 0
        assert ok,"parse of canonical smiles %s failed"%(smi)

        assert canonSmi==newCanon,'smiles canonicalization failed for %s\n%s != %s'%(smi,canonSmi,newCanon)

        pkl = cPickle.dumps(m)
        m2 = cPickle.loads(pkl)
        try3 = rdmol.MolToSmiles(m2)
        assert canonSmi == try3,'Non-canonical: mol %s gave %s (should be %s) after pickling'%(smi,try3,canonSmi)

        nFrags = max(rdmol.FindMolFrags(m))+1
        rings = rdmol.FindSSSR(m)
        cyclomat = m.getNumBonds() - m.getNumAtoms() + nFrags;
        assert cyclomat == len(rings),'bad num rings for %s\n\t%s!=%s'%(smi,cyclomat,len(rings))
        nChords = m.getNumBonds() - len(rdmol.FindSpanningTree(m))
        assert len(rings) == nChords, 'bad num chords for %s\n\t%s!=%s'%(smi,nChords,len(rings))
      
    
  def testBulkSmiles1(self):
    """ bulk processing of easy_smiles.txt """
    inLines = open('%s/easy_smiles.txt'%(self.dataPath)).readlines()
    self._bulkSmilesTest(inLines)

  def testBulkSmiles2(self):
    """ bulk processing of ntp_smiles.txt """
    inLines = open('%s/ntp_smiles.txt'%(self.dataPath)).readlines()
    self._bulkSmilesTest(inLines)

  def testBulkSmiles3(self):
    """ bulk processing of rtecs_smiles.txt """
    inLines = open('%s/rtecs_smiles.txt'%(self.dataPath)).readlines()
    self._bulkSmilesTest(inLines)

  def testSmilesProblems(self):
    """ some SMILES which have been problems at times """

    smi = 'Sc1ccc1'
    m = rdmol.MolFromSmiles(smi)
    assert m.getNumAtoms()==5,'bad numAtoms for: %s'%smi

    smi = '[Sc]1ccc1'
    m = rdmol.MolFromSmiles(smi)
    assert m.getNumAtoms()==4,'bad numAtoms for: %s'%smi

    smi = 'Sc'
    m = rdmol.MolFromSmiles(smi)
    assert m.getNumAtoms()==2,'bad numAtoms for: %s'%smi

    smi = 'C(=C)C'
    m = rdmol.MolFromSmiles(smi)
    assert m.getNumAtoms()==3,'bad numAtoms for: %s'%smi
    assert m.getBondWithIdx(0).getBondType() == rdmol.Bond.DOUBLE,'bad bond type for %s'%smi

    smi = 'C1=C(C1)'
    m = rdmol.MolFromSmiles(smi)
    assert m.getNumAtoms()==3,'bad numAtoms for: %s'%smi
    assert m.getBondWithIdx(0).getBondType() == rdmol.Bond.DOUBLE,'bad bond type for %s'%smi
    rings = rdmol.FindSSSR(m)
    assert len(rings)==1,'bad SSSR for %s'%smi
    assert len(rings[0])==3,'bad SSSR for %s'%smi

    smi = 'C(=C1)(C1)'
    m = rdmol.MolFromSmiles(smi)
    assert m.getNumAtoms()==3,'bad numAtoms for: %s'%smi
    assert m.getBondWithIdx(0).getBondType() == rdmol.Bond.DOUBLE,'bad bond type for %s'%smi
    rings = rdmol.FindSSSR(m)
    assert len(rings)==1,'bad SSSR for %s'%smi
    assert len(rings[0])==3,'bad SSSR for %s'%smi
    b = rdmol.GetBondBetweenAtoms(m,0,1)
    assert b,'gbba failed for %s'%smi
    assert b.getBondType()==rdmol.Bond.DOUBLE,'gbba order failed for %s'%smi
    b = rdmol.GetBondBetweenAtoms(m,1,2)
    assert b,'gbba failed for %s'%smi
    assert b.getBondType()==rdmol.Bond.SINGLE,'gbba order failed for %s'%smi
    b = rdmol.GetBondBetweenAtoms(m,2,0)
    assert b,'gbba failed for %s'%smi
    assert b.getBondType()==rdmol.Bond.SINGLE,'gbba order failed for %s'%smi
    b = rdmol.GetBondBetweenAtoms(m,1,0)
    assert b,'gbba failed for %s'%smi
    assert b.getBondType()==rdmol.Bond.DOUBLE,'gbba order failed for %s'%smi
    b = rdmol.GetBondBetweenAtoms(m,2,1)
    assert b,'gbba failed for %s'%smi
    assert b.getBondType()==rdmol.Bond.SINGLE,'gbba order failed for %s'%smi
    b = rdmol.GetBondBetweenAtoms(m,0,2)
    assert b,'gbba failed for %s'%smi
    assert b.getBondType()==rdmol.Bond.SINGLE,'gbba order failed for %s'%smi

  def testGetBonds(self):
    """ testing bond getters """
    smi = 'CCCCC'
    m = rdmol.MolFromSmiles(smi)
    for i in range(m.getNumAtoms()):
      for j in range(m.getNumAtoms()):
        b = rdmol.GetBondBetweenAtoms(m,i,j)
        if abs(i-j) == 1:
          assert b,'bad gbba'
          assert b.getBondType()==rdmol.Bond.SINGLE,"bad gbba order"
        else:
          assert not b,"improper gbba"

    
  def testFrags(self):
    """ test detecting fragments """

    smi = '[Na+].[Cl-]'
    m = rdmol.MolFromSmiles(smi)
    frags = list(rdmol.FindMolFrags(m))
    assert max(frags)==1,'bad numFrags for: %s'%smi
    assert frags.count(0) == 1,'bad frag0 count for: %s'%smi
    assert frags.count(1) == 1,'bad frag1 count for: %s'%smi

    smi = '[Na+].CC(=O)[O-]'
    m = rdmol.MolFromSmiles(smi)
    frags = list(rdmol.FindMolFrags(m))
    assert max(frags)==1,'bad numFrags for: %s'%smi
    assert frags.count(0) == 1,'bad frag0 count for: %s'%smi
    assert frags.count(1) == 4,'bad frag1 count for: %s'%smi

    smi = '[Na+].[Cl-].CCOCC.C1CCC1'
    m = rdmol.MolFromSmiles(smi)
    frags = list(rdmol.FindMolFrags(m))
    assert max(frags)==3,'bad numFrags for: %s'%smi
    assert frags.count(0) == 1,'bad frag0 count for: %s'%smi
    assert frags.count(1) == 1,'bad frag1 count for: %s'%smi
    assert frags.count(2) == 5,'bad frag2 count for: %s'%smi
    assert frags.count(3) == 4,'bad frag3 count for: %s'%smi


  def testSmi2CDXMLChains(self):
    """ testing some chains (tests SMILES parsing and CDXML generation) """
    d = [('CCCC','chain1.cdxml'),
         ('C=CC#C','chain2.cdxml'),
         ('COCC','chain3.cdxml'),
         ('CO>CC','chain4.cdxml'),
         ]
    for smi,cdxFile in d:
      cdxFile = '%s/%s'%(self.dataPath,cdxFile)
      refD = open(cdxFile,'r').read()
      genD = Smi2CDXML(smi)
      assert genD == refD, 'bad CDXML for SMILES: %s, file %s'%(smi,cdxFile)

  def testSmi2CDXMLRings(self):
    """ testing some rings (tests SMILES parsing and CDXML generation) """
    d = [('C1CCC1','ring1.cdxml'),
         ('C1CCC=1','ring2.cdxml'),
         ('c1ccccc1','ring3.cdxml'),
         ]
    for smi,cdxFile in d:
      cdxFile = '%s/%s'%(self.dataPath,cdxFile)
      refD = open(cdxFile,'r').read()
      genD = Smi2CDXML(smi)
      assert genD == refD, 'bad CDXML for SMILES: %s, file %s'%(smi,cdxFile)

  def testCDXMLParse1(self):
    """ tests cdxml parsing and various molecular getters """
    smi,cdxFile = ('CCCC','chain1.cdxml')
    cdxFile = '%s/%s'%(self.dataPath,cdxFile)
    cdxD = open(cdxFile,'r').read()
    m = rdmol.MolFromCDXML(cdxD)
    assert m.getNumAtoms()==4,"bad numAtoms"
    assert m.getNumBonds()==3,"bad numBonds"
    assert m.getBondWithIdx(0).getBondType()==Bond.SINGLE,"bad bond type"
    assert m.getBondWithIdx(1).getBondType()==Bond.SINGLE,"bad bond type"
    assert m.getBondWithIdx(2).getBondType()==Bond.SINGLE,"bad bond type"
    b = m.getBondWithIdx(2)
    assert b.getBeginAtomIdx()==2,"bad begin atom"
    assert b.getEndAtomIdx()==3,"bad end atom"
    assert b.getOtherAtomIdx(2)==3,"bad other atom"
    assert b.getOtherAtomIdx(3)==2,"bad other atom"

  def testCDXMLParse2(self):
    """ tests cdxml parsing and various molecular getters """
    smi,cdxFile = ('C=CC#C','chain2.cdxml')
    cdxFile = '%s/%s'%(self.dataPath,cdxFile)
    cdxD = open(cdxFile,'r').read()
    m = rdmol.MolFromCDXML(cdxD)
    assert m.getNumAtoms()==4,"bad numAtoms"
    assert m.getNumBonds()==3,"bad numBonds"
    assert m.getBondWithIdx(0).getBondType()==Bond.DOUBLE,"bad bond type"
    assert m.getBondWithIdx(1).getBondType()==Bond.SINGLE,"bad bond type"
    assert m.getBondWithIdx(2).getBondType()==Bond.TRIPLE,"bad bond type"
    b = m.getBondWithIdx(2)
    assert b.getBeginAtomIdx()==2,"bad begin atom"
    assert b.getEndAtomIdx()==3,"bad end atom"
    assert b.getOtherAtomIdx(2)==3,"bad other atom"
    assert b.getOtherAtomIdx(3)==2,"bad other atom"

  def testCDXMLParse3(self):
    """ tests cdxml parsing and various molecular getters """
    smi,cdxFile = ('C1CCC1','ring1.cdxml')
    cdxFile = '%s/%s'%(self.dataPath,cdxFile)
    cdxD = open(cdxFile,'r').read()
    m = rdmol.MolFromCDXML(cdxD)
    assert m.getNumAtoms()==4,"bad numAtoms"
    assert m.getNumBonds()==4,"bad numBonds"

    nbrs = rdmol.GetAtomNeighbors(m,m.getAtomWithIdx(0))
    assert len(nbrs) == 2,"bad len(nbrs)"
    assert 1 in nbrs,"bad nbrs"
    assert 3 in nbrs,"bad nbrs"
    for i in range(m.getNumAtoms()):
      assert len(rdmol.GetAtomNeighbors(m,m.getAtomWithIdx(i)))==2,"bad len(nbrs)"

    assert m.getBondWithIdx(0).getBondType()==Bond.SINGLE,"bad bond type"
    assert m.getBondWithIdx(1).getBondType()==Bond.SINGLE,"bad bond type"
    assert m.getBondWithIdx(2).getBondType()==Bond.SINGLE,"bad bond type"
    assert m.getBondWithIdx(3).getBondType()==Bond.SINGLE,"bad bond type"
    b = m.getBondWithIdx(2)
    assert b.getBeginAtomIdx()==2,"bad begin atom"
    assert b.getEndAtomIdx()==3,"bad end atom"
    assert b.getOtherAtomIdx(2)==3,"bad other atom"
    assert b.getOtherAtomIdx(3)==2,"bad other atom"
    
    rings = rdmol.FindSSSR(m)
    assert len(rings)==1,"bad SSSR count"
    assert len(rings[0])==4,"bad SSSR ring size"
    for idx in range(m.getNumAtoms()):
      assert idx in rings[0],"atom not in SSSR"

    # test the property caching stuff
    p = rdmol.MolProps(m)
    rings = p.FindSSSR()
    assert len(rings)==1,"bad SSSR count"
    assert len(rings[0])==4,"bad SSSR ring size"
    for idx in range(m.getNumAtoms()):
      assert idx in rings[0],"atom not in SSSR"



  def testCDXMLParse4(self):
    """ tests cdxml parsing and various molecular getters """
    smi,cdxFile = ('c1ccccc11','ring3.cdxml')
    cdxFile = '%s/%s'%(self.dataPath,cdxFile)
    cdxD = open(cdxFile,'r').read()
    m = rdmol.MolFromCDXML(cdxD)
    assert m.getNumAtoms()==6,"bad numAtoms"
    assert m.getNumBonds()==6,"bad numBonds"

    nbrs = rdmol.GetAtomNeighbors(m,m.getAtomWithIdx(0))
    assert len(nbrs) == 2,"bad len(nbrs)"
    assert 1 in nbrs,"bad nbrs"
    assert 2 not in nbrs,"bad nbrs"
    assert 5 in nbrs,"bad nbrs"
    for i in range(m.getNumAtoms()):
      assert len(rdmol.GetAtomNeighbors(m,m.getAtomWithIdx(i)))==2,"bad len(nbrs)"

    for i in range(m.getNumBonds()):
      assert m.getBondWithIdx(i).getBondType()==Bond.AROMATIC,"bad bond type"
    b = m.getBondWithIdx(2)
    assert b.getBeginAtomIdx()==2,"bad begin atom"
    assert b.getEndAtomIdx()==3,"bad end atom"
    assert b.getOtherAtomIdx(2)==3,"bad other atom"
    assert b.getOtherAtomIdx(3)==2,"bad other atom"
    
    rings = rdmol.FindSSSR(m)
    assert len(rings)==1,"bad SSSR count"
    assert len(rings[0])==6,"bad SSSR ring size"
    for idx in range(m.getNumAtoms()):
      assert idx in rings[0],"atom not in SSSR"

    # test the property caching stuff
    p = rdmol.MolProps(m)
    rings = p.FindSSSR()
    assert len(rings)==1,"bad SSSR count"
    assert len(rings[0])==6,"bad SSSR ring size"
    for idx in range(m.getNumAtoms()):
      assert idx in rings[0],"atom not in SSSR"


  def testAtomOps(self):
    """ testing atoms """
    m = Mol()
    nums = [6,8,7]
    symbs = ['C','O','N']
    a = Atom(nums[0])
    assert a.getAtomicNum()==nums[0],"bad atomic number"

    m.addAtom(a)
    assert m.getNumAtoms()==1,"bad num atoms" 
    m.addAtom(Atom(nums[1]))
    assert m.getNumAtoms()==2,"bad num atoms" 
    m.addAtom(Atom(nums[2]))
    for i in range(m.getNumAtoms()):
      a = m.getAtomWithIdx(i)
      assert a.getIdx() == i,"bad atom idx"
      assert a.getSymbol() == symbs[i],'bad symbol'
      assert a.getAtomicNum() == nums[i],"bad atomic num"
      assert a.getDegree() == 0, "bad degree"
      assert a.getExplicitValence() == 0, "bad e valence: %d"%a.getExplicitValence()
      assert a.getFormalCharge() == 0, "bad chg"

      a.setNoImplicit(1)
      a = m.getAtomWithIdx(i)
      assert a.getNoImplicit() == 1, "bad noimplicit"
      assert a.getImplicitValence() == 0, "bad i valence"
      assert a.getNumExplicitHs() == 0, "bad # e Hs"
      assert a.getIsAromatic() == 0, "bad aromatic"
      

  def testSubstructPass(self):
    """ testing substructure matches which should pass """
    data = [('C1CCCCC1','CCC'),('C1CCCCC1','CCCCC'),
            ('CCC(CC)CC','CC(C)C')]
    for molSmi,querySmi in data:
      mol = rdmol.MolFromSmiles(molSmi)
      query = rdmol.MolFromSmiles(querySmi)
      assert rdmol.HasSubstructMatch(mol,mol),"no mol self match"
      assert rdmol.HasSubstructMatch(mol,query),"no match"
      assert rdmol.HasSubstructMatch(query,query),"no query self match"


  def testAdjMat(self):
    smi = 'CC'
    mol = rdmol.MolFromSmiles(smi)
    tgt = array([[0,1],[1,0]])
    mat = rdmol.GetAdjMat(mol,0)
    assert _fmatcmp(mat,tgt) < _FTOL,'bad adj mat for %s'%(smi)
    tgt = array([[0,1],[1,0]])
    mat = rdmol.GetAdjMat(mol,1)
    assert _fmatcmp(mat,tgt) < _FTOL,'bad val adj mat for %s'%(smi)
    
    smi = 'C=C'
    mol = rdmol.MolFromSmiles(smi)
    tgt = array([[0,1],[1,0]])
    mat = rdmol.GetAdjMat(mol,0)
    assert _fmatcmp(mat,tgt) < _FTOL,'bad adj mat for %s'%(smi)
    tgt = array([[0,2],[2,0]])
    mat = rdmol.GetAdjMat(mol,1)
    assert _fmatcmp(mat,tgt) < _FTOL,'bad val adj mat for %s'%(smi)
    
    smi = 'C#C'
    mol = rdmol.MolFromSmiles(smi)
    tgt = array([[0,1],[1,0]])
    mat = rdmol.GetAdjMat(mol,0)
    assert _fmatcmp(mat,tgt) < _FTOL,'bad adj mat for %s'%(smi)
    tgt = array([[0,3],[3,0]])
    mat = rdmol.GetAdjMat(mol,1)
    assert _fmatcmp(mat,tgt) < _FTOL,'bad val adj mat for %s'%(smi)
    
    smi = 'C.C'
    mol = rdmol.MolFromSmiles(smi)
    tgt = array([[0,0],[0,0]])
    mat = rdmol.GetAdjMat(mol,0)
    assert _fmatcmp(mat,tgt) < _FTOL,'bad adj mat for %s'%(smi)
    tgt = array([[0,0],[0,0]])
    mat = rdmol.GetAdjMat(mol,1)
    assert _fmatcmp(mat,tgt) < _FTOL,'bad val adj mat for %s'%(smi)
    
    smi = 'cc'
    mol = rdmol.MolFromSmiles(smi)
    tgt = array([[0,1],[1,0]])
    mat = rdmol.GetAdjMat(mol,0)
    assert _fmatcmp(mat,tgt) < _FTOL,'bad adj mat for %s'%(smi)
    tgt = array([[0,1.5],[1.5,0]])
    mat = rdmol.GetAdjMat(mol,1)
    assert _fmatcmp(mat,tgt) < _FTOL,'bad val adj mat for %s'%(smi)
    
    smi = 'C>C'
    mol = rdmol.MolFromSmiles(smi)
    tgt = array([[0,1],[1,0]])
    mat = rdmol.GetAdjMat(mol,0)
    assert _fmatcmp(mat,tgt) < _FTOL,'bad adj mat for %s'%(smi)
    tgt = array([[0,0],[1,0]])
    mat = rdmol.GetAdjMat(mol,1)
    assert _fmatcmp(mat,tgt) < _FTOL,'bad val adj mat for %s: %s'%(smi,str(mat))
    
    smi = 'C<C'
    mol = rdmol.MolFromSmiles(smi)
    tgt = array([[0,1],[1,0]])
    mat = rdmol.GetAdjMat(mol,0)
    assert _fmatcmp(mat,tgt) < _FTOL,'bad adj mat for %s'%(smi)
    tgt = array([[0,1],[0,0]])
    mat = rdmol.GetAdjMat(mol,1)
    assert _fmatcmp(mat,tgt) < _FTOL,'bad val adj mat for %s: %s'%(smi,str(mat))
    

      
  def testAtomListPass(self):
    """ testing atom list matches which should pass """
    smis = ['C1=CC=CC=C1','C1C=CC=CC=1',
            'N1=CC=CC=C1','C1=NC=CC=C1','C1=CN=CC=C1',
            'C1=CC=NC=C1','C1=CC=CN=C1','C1=CC=CC=N1',
            'P1=CC=CC=C1','C1=PC=CC=C1','C1=CP=CC=C1',
            'C1=CC=PC=C1','C1=CC=CP=C1','C1=CC=CC=P1',
            'C1=C(C)C=CC=C1','C1C=C(CC)C=C(C)C=1',
            ]
    cdxFile = '%s/list-query.cdxml'%(self.dataPath)
    query = rdmol.MolFromCDXML(open(cdxFile,'r').read())
    for smi in smis:
      mol = rdmol.MolFromSmiles(smi)
      assert rdmol.HasSubstructMatch(mol,query),"no match found"
      matches = rdmol.SubstructMatch(mol,query)
      assert len(matches[0]) == 6, "bad match length"

  def testAtomListFail(self):
    """ testing atom list matches which should fail """
    smis = ['C1CC=CC=C1','c1ccccc1','c1ccccc1','c1ccccc1','c1ccccc1',
            'C1=CP=CO=C1','O1=CC=CC=C1','C1=NC=CN=C1',
            ]
    cdxFile = '%s/list-query.cdxml'%(self.dataPath)
    query = rdmol.MolFromCDXML(open(cdxFile,'r').read())
    for smi in smis:
      mol = rdmol.MolFromSmiles(smi)
      assert not rdmol.HasSubstructMatch(mol,query),"invalid match found"

  def testBondListPass(self):
    """ testing bond list matches which should pass """
    smis = ['CCCC=C','C=CCCC','CC=CC=C','C=CC=CC',
            'C1CCC=CC1','C1CC=CC=C1'
            'C=CC=CCCCO','CC=C(C=C)COC',
            ]
    cdxFile = '%s/bond-query.cdxml'%(self.dataPath)
    query = rdmol.MolFromCDXML(open(cdxFile,'r').read())
    for smi in smis:
      mol = rdmol.MolFromSmiles(smi)
      assert rdmol.HasSubstructMatch(mol,query),"no match found"
      matches = rdmol.SubstructMatch(mol,query)
      assert len(matches[0]) == 5, "bad match length"

  def testBondListFail(self):
    """ testing bond list matches which should fail """
    smis = ['CCCCC','C=COCC',
            ]
    cdxFile = '%s/bond-query.cdxml'%(self.dataPath)
    query = rdmol.MolFromCDXML(open(cdxFile,'r').read())
    for smi in smis:
      mol = rdmol.MolFromSmiles(smi)
      assert not rdmol.HasSubstructMatch(mol,query),"invalid match found"

  def testMolParse(self):
    """ test basic mol file parsing """
    molFile = '%s/mol1.mol'%(self.dataPath)
    m = rdmol.MolFromMolFile(molFile)

    assert m.getNumAtoms() == 5,'bad numAtoms'
    assert rdmol.MolToSmiles(m) == '[C-1]1C=CC=C1'

    mP = rdmol.MolProps(m)
    sssr = mP.FindSSSR()
    assert len(sssr) == 1, 'bad sssr'
    assert len(sssr[0]) == 5, 'bad sssr'

  def testMolAtomListPass(self):
    """ testing atom queries from Mol files """
    smis = ['C1=CC=CC=C1','C1C=CC=CC=1',
        'N1=CC=CC=C1','C1=NC=CC=C1','C1=CN=CC=C1',
        'C1=CC=NC=C1','C1=CC=CN=C1','C1=CC=CC=N1',
        'P1=CC=CC=C1','C1=PC=CC=C1','C1=CP=CC=C1',
        'C1=CC=PC=C1','C1=CC=CP=C1','C1=CC=CC=P1',
        'C1=C(C)C=CC=C1','C1C=C(CC)C=C(C)C=1',
        ]
    molFile = '%s/list-query.mol'%(self.dataPath)
    query = rdmol.MolFromMolFile(molFile)
    for smi in smis:
      mol = rdmol.MolFromSmiles(smi)
      assert rdmol.HasSubstructMatch(mol,query),"match not found"

    molFile = '%s/query2.mol'%(self.dataPath)
    query = rdmol.MolFromMolFile(molFile)
    for smi in smis:
      mol = rdmol.MolFromSmiles(smi)
      assert rdmol.HasSubstructMatch(mol,query),"match not found"

  def testMolAtomListFail(self):
    """ testing atom list matches which should fail """
    smis = ['C1CC=CC=C1','c1ccccc1','c1ccccc1','c1ccccc1','c1ccccc1',
            'C1=CP=CO=C1','O1=CC=CC=C1','C1=NC=CN=C1',
            ]
    molFile = '%s/list-query.mol'%(self.dataPath)
    query = rdmol.MolFromMolFile(molFile)
    for smi in smis:
      mol = rdmol.MolFromSmiles(smi)
      assert not rdmol.HasSubstructMatch(mol,query),"invalid match found"

    molFile = '%s/query2.mol'%(self.dataPath)
    query = rdmol.MolFromMolFile(molFile)
    for smi in smis:
      mol = rdmol.MolFromSmiles(smi)
      assert not rdmol.HasSubstructMatch(mol,query),"invalid match found"


  def testMolBondListPass(self):
    """ testing bond list matches which should pass """
    # query is CC?CC=C where ? is single or double
    smis = ['CCCC=C','C=CCCC','CC=CC=C','C=CC=CC',
            'C1CCC=CC1','C1CC=CC=C1'
            'C=CC=CCCCO','CC=C(C=C)COC',
            ]
    molFile = '%s/bond-query.mol'%(self.dataPath)
    query = rdmol.MolFromMolFile(molFile)
    for smi in smis:
      mol = rdmol.MolFromSmiles(smi)
      assert rdmol.HasSubstructMatch(mol,query),"no match found"
      matches = rdmol.SubstructMatch(mol,query)
      assert len(matches[0]) == 5, "bad match length"

  def testMolBondListFail(self):
    """ testing bond list matches which should fail """
    # query is CC?CC=C where ? is single or double
    smis = ['CCCCC','C=COCC','CC#CC=C','CCCC#C',
            ]
    molFile = '%s/bond-query.mol'%(self.dataPath)
    query = rdmol.MolFromMolFile(molFile)
    for smi in smis:
      mol = rdmol.MolFromSmiles(smi)
      assert not rdmol.HasSubstructMatch(mol,query),"invalid match found"
    
  def testMolBondList2(self):
    """ testing additional bond list matches """
    # query is CC?CC=C where ? is single or aromatic
    smis = ['CCCC=C','C=CCCC','c1ccc(C)c(C=C)c1',
            ]
    molFile = '%s/bond-query2.mol'%(self.dataPath)
    query = rdmol.MolFromMolFile(molFile)
    for smi in smis:
      mol = rdmol.MolFromSmiles(smi)
      assert rdmol.HasSubstructMatch(mol,query),"no match found"
      matches = rdmol.SubstructMatch(mol,query)
      assert len(matches[0]) == 5, "bad match length"

  def testMolBondList3(self):
    """ testing additional bond list matches """
    # query is CC?CC=C where ? is double or aromatic
    smis = ['CC=CC=C','C=CC=CC','c1ccc(C)c(C=C)c1',
            ]
    molFile = '%s/bond-query3.mol'%(self.dataPath)
    query = rdmol.MolFromMolFile(molFile)
    for smi in smis:
      mol = rdmol.MolFromSmiles(smi)
      assert rdmol.HasSubstructMatch(mol,query),"no match found"
      matches = rdmol.SubstructMatch(mol,query)
      assert len(matches[0]) == 5, "bad match length"

  def testMolBondList4(self):
    """ testing additional bond list matches """
    # query is CC?CC=C where ? is any type of bond
    smis = ['CC=CC=C','C=CC=CC','c1ccc(C)c(C=C)c1','C=CC#CC','C=CCCC',
            ]
    molFile = '%s/bond-query4.mol'%(self.dataPath)
    query = rdmol.MolFromMolFile(molFile)
    for smi in smis:
      mol = rdmol.MolFromSmiles(smi)
      assert rdmol.HasSubstructMatch(mol,query),"no match found"
      matches = rdmol.SubstructMatch(mol,query)
      assert len(matches[0]) == 5, "bad match length"

  def _smiTest(self,smiList):
    for smi,spellings in smiList:
      ok = 1
      try:
        m = rdmol.MolFromSmiles(smi)
      except:
        ok = 0
      assert ok,'parse failure for SMILES: %s'%smi  
      canSmi = rdmol.MolToSmiles(m)
      for spelling in spellings:
        m = rdmol.MolFromSmiles(spelling)
        trySmi = rdmol.MolToSmiles(m)
        assert canSmi==trySmi,'Non-canonical: mol %s gave %s (should be %s)'%(spelling,trySmi,canSmi)
        m2 = rdmol.MolFromSmiles(trySmi)
        try2 = rdmol.MolToSmiles(m2)
        assert canSmi==try2,'Non-canonical: mol %s gave %s (should be %s) on second pass'%(spelling,try2,canSmi)
      pkl = cPickle.dumps(m)
      m2 = cPickle.loads(pkl)
      try3 = rdmol.MolToSmiles(m2)
      assert canSmi == try3,'Non-canonical: mol %s gave %s (should be %s) after pickling'%(smi,try3,canSmi)
        
    
  def testLinear1(self):
    " testing first batch of linear mols "
    smiList = [('O=CCO',('OCC=O',
                         'C(O)C=O',
                         'C(C=O)O',
                         'C(CO)=O',
                         )),
               ('OCC(C=C)CCC(C#N)CC',('C=CC(CO)CCC(C#N)CC',
                                      'C(CO)(C=C)CCC(CC)C#N',
                                      'C(CO)(C=C)CCC(C#N)CC',
                                      'C(C=C)(CO)CCC(C#N)CC',
                                      'C(C=C)(CO)CCC(CC)C#N',
                                      )),
               ('[Se]=CCO',('OCC=[Se]',
                         'C(O)C=[Se]',
                         'C(C=[Se])O',
                         'C(CO)=[Se]',
                         )),
               ('O=C(OC)C',('COC(=O)C','O=C(C)OC','O(C)C(=O)C')),
               ]
    self._smiTest(smiList)
  def testRings1(self):
    " testing first batch of rings "
    smiList = [('C1OCCCC1',('O1CCCCC1',
                            'C1COCCC1',
                            'C1CCOCC1',
                            'C1CCCOC1',
                            'C1CCCCO1',)),
               ('CC1=CCCCC1',('C1=C(C)CCCC1',
                              'C1CC=C(C)CC1',)),
               ('CC1C=CCCC1',('C1=CC(C)CCC1',
                              'C1CC=CC(C)C1',)),
               ('C1=CC=CC=C1',('C=1C=CC=CC1',
                               'C1C=CC=CC=1')),
               ]
    self._smiTest(smiList)
        
  def testRings2(self):
    " testing second batch of rings "
    smiList = [('c1c(cc2nc3cc(ccc3cc2c1))',('c1ccc2cc3ccccc3nc2c1',
                                            'c1ccc2nc3ccccc3cc2c1',
                                            'c1c2nc3ccccc3cc2ccc1')),
               ('Cc1ccc2nc3ccccc3cc2c1',('c1ccc2nc3ccc(C)cc3cc2c1',
                                          )),
               ('c1c(C)cc2nc3ccccc3cc2c1',('c1ccc2nc3cc(C)ccc3cc2c1',
                                           )),

               ]
    self._smiTest(smiList)

  def testProblems(self):
    " testing molecules which have been problematic "

    smiList = [ ('[Al+3]CCC',('CCC[Al+3]','C(C)(C[Al+3])' ) ),
                ('C(=O)(Cl)CC(=O)Cl',('ClC(CC(Cl)=O)=O','C(Cl)(=O)CC(=O)Cl','C(Cl)(=O)CC(Cl)=O')),
#                ('C(=O)(Cl)c1ccc(C(=O)Cl)cc1',('O=C(Cl)c1ccc(cc1)C(Cl)=O','C(Cl)(=O)C1=CC=C(C=C1)C(Cl)=O','ClC(=O)c1ccc(cc1)C(=O)Cl')),
#                ('[N+](=O)([O-])c1ccc([N+](=O)[O-])cc1',('[N+]([O-])(=O)C1=CC=C(C=C1)[N+](=O)[O-]','O=[N+1]([O--1])c1ccc(cc1)[N+1]([O--1])=O','[O--1][N+1](=O)c1ccc(cc1)[N+1]([O--1])=O')),
                ('Oc1c3c(cc(c1)S(=O)(=O)O)cc(NC(=O)c2ccccc2)cc3',
                  ('O=S(=O)(O)c1cc(O)c2ccc(NC(=O)c3ccccc3)cc2c1',
                  'OS(=O)(=O)c1cc(O)c2ccc(NC(=O)c3ccccc3)cc2c1')),
                ('C(C(C)(C)O)C(C)O',('C([C@@](C)(O)C)C(C)O','CC(O)CC(C)(O)C','CC(O)CC(O)(C)C')),
                ('C',('C')),
                ('C(Cl)(Br)(F)CC(Cl)(Br)(F)',('C(Cl)(F)(Br)CC(F)(Cl)(Br)','C(Cl)(Br)(F)CC(Cl)(F)(Br)',
                                              'C(F)(Br)(Cl)CC(Br)(Cl)(F)','C(C(Cl)(Br)(F))C(F)(Cl)Br')),
                  
               ]
    self._smiTest(smiList)
        
  def testHighSymmetry(self):
    " testing tricky (high-symmetry) molecules "

    smiList = [('CC(C)CC',('CCC(C)C',)),
               ('C1CCCC1CCC',('CCCC1CCCC1',)),
               ('C1(C)CC(C)CCC1',('CC1CCCC(C)C1',)),
               ('CCC1CCCCC1CC',('CCC1CCCCC1CC',)),
               ('CCC1CC(CC)CCC1',('CCC1CCCC(CC)C1',)),
               ('C1CCCCC1CC(CC)CC',('CCC(CC)CC1CCCCC1',)),
               ('C1CCCC2C1CC(CC)CC2',('CCC1CCC2CCCCC2C1',)),
               ('CC1CCCC2C1C(C)CCC2',('CC1CCCC2CCCC(C)C12',)),
               ('C2CCC1CCC(C)C12',('CC1CCC2CCCC12',)),
               ('CC(C)CCCC(C)C',('CC(CCCC(C)C)C',)),
               ]
    self._smiTest(smiList)

  def testFailures(self):
    " testing molecules which are known to fail (EXPECT FAILURES)"

    smiList = [('C13C6C1C2C4C2C3C5C4C56',
                ('C45C1C6C3C6C5C4C2C3C12',
                 'C45C2C6C3C6C5C4C1C3C12')),
                ]
    self._smiTest(smiList)


  def testDblBonds(self):
    " testing double bonds "
    smiList = [(r'C/C=C/O',(r'C\C=C\O',r'O/C=C/C',r'O\C=C\C')),
               (r'C/C=C\O',(r'C\C=C/O',r'O\C=C/C',r'O/C=C\C')),
               (r'C/C=C\O',(r'C\C=C/O',r'O\C=C/C',r'O/C=C\C')),
               (r'C/C(O)=C\O',(r'C\C(O)=C/O',r'O\C=C(O)/C',r'O/C=C(\C)O')),
               (r'C/C(O)=C(\O)C',(r'O\C(C)=C(\O)C',r'O\C(C)=C(/C)O',)),
               ]
    self._smiTest(smiList)
    

  def testDativeBonds(self):
    " testing dative bonds "
    smiList = [(r'S<CO',('S<(CO)','S(<CO)','OC>S','C(O)>S')),
               ]
    self._smiTest(smiList)
    
    m = rdmol.MolFromSmiles('S<CO')
    assert m.getAtomWithIdx(0).getExplicitValence()==1,'bad explicit valence'
    assert m.getAtomWithIdx(1).getExplicitValence()==1,'bad explicit valence'
    assert m.getAtomWithIdx(2).getExplicitValence()==1,'bad explicit valence'

  def _smiTest2(self,smiList):
    # don't test re-parsing smiles since we don't deal with writing hapto things
    #  (yet)
    for smi,spellings in smiList:
      ok = 1
      try:
        m = rdmol.MolFromSmiles(smi)
      except:
        ok = 0
      assert ok,'parse failure for SMILES: %s'%smi  
      canSmi = rdmol.MolToSmiles(m)
      for spelling in spellings:
        m = rdmol.MolFromSmiles(spelling)
        trySmi = rdmol.MolToSmiles(m)
        assert canSmi==trySmi,'Non-canonical: mol %s gave %s (should be %s)'%(spelling,trySmi,canSmi)
      pkl = cPickle.dumps(m)
      m2 = cPickle.loads(pkl)
      try3 = rdmol.MolToSmiles(m2)
      assert canSmi == try3,'Non-canonical: mol %s gave %s (should be %s) after pickling'%(smi,try3,canSmi)

  def testHaptoParse(self):
    " testing parsing of hapto bonds "

    # different writings of the same thing
    smi = r'[Fe](<_(c_1c_c_c_c_c_1))(C)C'
    m = rdmol.MolFromSmiles(smi)
    fe = m.getAtomWithIdx(0)
    c1 = m.getAtomWithIdx(1)
    c2 = m.getAtomWithIdx(6)
    assert fe.getExplicitValence()==8,'bad explicit valence for %s'%smi
    assert c1.getExplicitValence()==3,'bad explicit valence for %s'%smi
    assert c2.getExplicitValence()==3,'bad explicit valence for %s'%smi
    smi = r'[Fe](C)(<_(c_1c_c_c_c_c_1))C'
    m = rdmol.MolFromSmiles(smi)
    fe = m.getAtomWithIdx(0)
    c1 = m.getAtomWithIdx(2)
    c2 = m.getAtomWithIdx(7)
    assert fe.getExplicitValence()==8,'bad explicit valence for %s'%smi
    assert c1.getExplicitValence()==3,'bad explicit valence for %s'%smi
    assert c2.getExplicitValence()==3,'bad explicit valence for %s'%smi
    smi = r'[Fe](C)(C)<_(c_1c_c_c_c_c_1)'
    m = rdmol.MolFromSmiles(smi)
    fe = m.getAtomWithIdx(0)
    c1 = m.getAtomWithIdx(3)
    c2 = m.getAtomWithIdx(8)
    assert fe.getExplicitValence()==8,'bad explicit valence for %s'%smi
    assert c1.getExplicitValence()==3,'bad explicit valence for %s'%smi
    assert c2.getExplicitValence()==3,'bad explicit valence for %s'%smi
    
    # a "slipped" ring
    smi = r'[Fe](<_(c_1c_c_ccc1))(C)C'
    m = rdmol.MolFromSmiles(smi)
    fe = m.getAtomWithIdx(0)
    c1 = m.getAtomWithIdx(1)
    c2 = m.getAtomWithIdx(6)
    assert fe.getExplicitValence()==5,'bad explicit valence for %s'%smi
    assert c1.getExplicitValence()==3,'bad explicit valence for %s'%smi
    assert c2.getExplicitValence()==3,'bad explicit valence for %s'%smi
    
    # two rings
    smi = r'[Fe](<_(c_1c_c_c_c_c_1))(<_(c_1c_c_c_c_c_1))(C)C'
    m = rdmol.MolFromSmiles(smi)
    fe = m.getAtomWithIdx(0)
    c1 = m.getAtomWithIdx(1)
    c2 = m.getAtomWithIdx(6)
    assert fe.getExplicitValence()==14,'bad explicit valence for %s'%smi
    assert c1.getExplicitValence()==3,'bad explicit valence for %s'%smi
    assert c2.getExplicitValence()==3,'bad explicit valence for %s'%smi

    # sandwich compound 1
    smi = r'[Fe](<_(c_1c_c_c_c_c_1))<_(c_1c_c_c_c_c_1)_>[Ti](Cl)(Cl)<_(c_1c_c_c_c_c_1)'
    m = rdmol.MolFromSmiles(smi)
    fe = m.getAtomWithIdx(0)
    ti = m.getAtomWithIdx(13)
    assert fe.getExplicitValence()==12,'bad explicit valence for %s'%smi
    assert ti.getExplicitValence()==14,'bad explicit valence for %s'%smi
    for i in range(m.getNumAtoms()):
      atom = m.getAtomWithIdx(i)
      if atom.getSymbol() == 'C':
        assert atom.getExplicitValence()==3,'bad C valence for %s'%smi

    

  def testHaptoCanon(self):
    " testing hapto canonicalization "
    smis = [(r'[Fe](<_(c_1c_c_c_c_c_1))(C)C',(r'[Fe](C)(<_(c_1c_c_c_c_c_1))C',
                                              r'[Fe](C)(C)<_(c_1c_c_c_c_c_1)',)),
            (r'[Fe](<_(c_1c_c_c_c_c_1))(<_(c_1c_c_c_c_c_1))(C)C',
             (r'[Fe](<_(c_1c_c_c_c_c_1))(C)(C)(<_(c_1c_c_c_c_c_1))',
              r'[Fe](C)(C)(<_(c_1c_c_c_c_c_1))(<_(c_1c_c_c_c_c_1))',
              r'[Fe](<_(c_1c_c_c_c_c_1))(C)(C)<_(c_1c_c_c_c_c_1)',
              r'[Fe](<_(c_1c_c_c_c_c_1))(C)(<_(c_1c_c_c_c_c_1))C',
              r'C[Fe](<_(c_1c_c_c_c_c_1))(C)(<_(c_1c_c_c_c_c_1))',
              r'C[Fe](<_(c_1c_c_c_c_c_1))(C)<_(c_1c_c_c_c_c_1)',
              ),),

            ]
    self._smiTest2(smis)
    

  def testMolPropsFrags(self):
    """ test detecting fragments """

    smi = '[Na+].[Cl-]'
    m = rdmol.MolFromSmiles(smi)
    p = rdmol.MolProps(m)
    frags = list(p.FindMolFrags())
    assert max(frags)==1,'bad numFrags for: %s'%smi
    assert frags.count(0) == 1,'bad frag0 count for: %s'%smi
    assert frags.count(1) == 1,'bad frag1 count for: %s'%smi
    # add an atom, this shouldn't change the props, because we're cached
    p.getMol().addAtom(rdmol.Atom(5))
    frags = list(p.FindMolFrags())
    assert max(frags)==1,'bad numFrags for: %s'%smi
    assert frags.count(0) == 1,'bad frag0 count for: %s'%smi
    assert frags.count(1) == 1,'bad frag1 count for: %s'%smi
    # reset the cache, so we recalculate
    p.reset()
    frags = list(p.FindMolFrags())
    assert max(frags)==2,'bad numFrags for: %s'%smi
    assert frags.count(0) == 1,'bad frag0 count for: %s'%smi
    assert frags.count(1) == 1,'bad frag1 count for: %s'%smi
    assert frags.count(2) == 1,'bad frag2 count for: %s'%smi


    
    smi = '[Na+].CC(=O)[O-]'
    m = rdmol.MolFromSmiles(smi)
    p = rdmol.MolProps(m)
    frags = list(p.FindMolFrags())
    assert max(frags)==1,'bad numFrags for: %s'%smi
    assert frags.count(0) == 1,'bad frag0 count for: %s'%smi
    assert frags.count(1) == 4,'bad frag1 count for: %s'%smi
    p.getMol().addAtom(rdmol.Atom(5))
    frags = list(p.FindMolFrags())
    assert max(frags)==1,'bad numFrags for: %s'%smi
    assert frags.count(0) == 1,'bad frag0 count for: %s'%smi
    assert frags.count(1) == 4,'bad frag1 count for: %s'%smi
    p.reset()
    frags = list(p.FindMolFrags())
    assert max(frags)==2,'bad numFrags for: %s'%smi
    assert frags.count(0) == 1,'bad frag0 count for: %s'%smi
    assert frags.count(1) == 4,'bad frag1 count for: %s'%smi
    assert frags.count(2) == 1,'bad frag2 count for: %s'%smi


    smi = '[Na+].[Cl-].CCOCC.C1CCC1'
    m = rdmol.MolFromSmiles(smi)
    p = rdmol.MolProps(m)
    frags = list(p.FindMolFrags())
    assert max(frags)==3,'bad numFrags for: %s'%smi
    assert frags.count(0) == 1,'bad frag0 count for: %s'%smi
    assert frags.count(1) == 1,'bad frag1 count for: %s'%smi
    assert frags.count(2) == 5,'bad frag2 count for: %s'%smi
    assert frags.count(3) == 4,'bad frag3 count for: %s'%smi
    p.getMol().addAtom(rdmol.Atom(5))
    frags = list(p.FindMolFrags())
    assert max(frags)==3,'bad numFrags for: %s'%smi
    assert frags.count(0) == 1,'bad frag0 count for: %s'%smi
    assert frags.count(1) == 1,'bad frag1 count for: %s'%smi
    assert frags.count(2) == 5,'bad frag2 count for: %s'%smi
    assert frags.count(3) == 4,'bad frag3 count for: %s'%smi
    p.reset()
    frags = list(p.FindMolFrags())
    assert max(frags)==4,'bad numFrags for: %s'%smi
    assert frags.count(0) == 1,'bad frag0 count for: %s'%smi
    assert frags.count(1) == 1,'bad frag1 count for: %s'%smi
    assert frags.count(2) == 5,'bad frag2 count for: %s'%smi
    assert frags.count(3) == 4,'bad frag3 count for: %s'%smi
    assert frags.count(4) == 1,'bad frag4 count for: %s'%smi
    



    
if __name__ == '__main__':
  unittest.main()

