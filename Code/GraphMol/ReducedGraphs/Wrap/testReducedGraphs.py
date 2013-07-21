# $Id$
# 
from rdkit import Chem
from rdkit.Chem import rdReducedGraphs as rdRG
from rdkit import RDConfig
import numpy
import unittest

class TestCase(unittest.TestCase) :
  def setUp(self):
    pass
  def test1(self):
    m = Chem.MolFromSmiles('OCCc1ccccc1')
    mrg = rdRG.GenerateMolExtendedReducedGraph(m)
    mrg.UpdatePropertyCache(False)
    self.failUnlessEqual('[*]cCCO',Chem.MolToSmiles(mrg))
    
    m = Chem.MolFromSmiles('OCCC1CCCCC1')
    mrg = rdRG.GenerateMolExtendedReducedGraph(m)
    mrg.UpdatePropertyCache(False)
    self.failUnlessEqual('[*]CCCO',Chem.MolToSmiles(mrg))
    
  def test2(self):
    m = Chem.MolFromSmiles('OCCc1ccccc1')
    mrg = rdRG.GenerateMolExtendedReducedGraph(m)
    mrg.UpdatePropertyCache(False)
    self.failUnlessEqual('[*]cCCO',Chem.MolToSmiles(mrg))

    fp1 = rdRG.GenerateErGFingerprintForReducedGraph(mrg)
    fp2 = rdRG.GetErGFingerprint(m)
    md = max(abs(fp1-fp2))
    self.failUnless(md<1e-4)


if __name__ == '__main__':
  unittest.main()
