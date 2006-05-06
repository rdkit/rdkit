import RDConfig
import unittest,os,sys,re
from Chem import rdchem_new as rdmol
from Chem.rdchem_new import Atom,Bond,Mol
from Numeric import *
import cPickle

class TestCase(unittest.TestCase):
  def setUp(self):
    print '\n%s: '%self.shortDescription()
    self.utilsPath = '%s/Code/GraphMol/Utils'%(RDConfig.RDBaseDir)
    self.substructPath = '%s/Code/GraphMol/Substruct'%(RDConfig.RDBaseDir)
    self.dataPath = '%s/Code/GraphMol/Wrap/test_data'%(RDConfig.RDBaseDir)
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
        
  def _bulkSmilesTest(self,lines):
    splitExpr = re.compile('[\t\ ]')
    for line in lines:
      if line[0] != '#':
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

    
  def testBulkSmiles3(self):
    """ bulk processing of rtecs_smiles.txt """
    inLines = open('%s/rtecs_smiles.5000.txt'%(self.dataPath)).readlines()
    self._bulkSmilesTest(inLines)


if __name__ == '__main__':
  unittest.main()

  






