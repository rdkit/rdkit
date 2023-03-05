import unittest

from rdkit import Chem
from rdkit.Chem import rdMolInterchange
from rdkit import RDConfig


class TestCase(unittest.TestCase):

  def setUp(self):
    pass

  def test1(self):
    smis = ('c1ccccc1','C[C@H](F)Cl')
    for smi in smis:
      m = Chem.MolFromSmiles(smi)
      csmi = Chem.MolToSmiles(m)
      json = rdMolInterchange.MolToJSON(m)
      nms = rdMolInterchange.JSONToMols(json)
      self.assertEqual(len(nms),1)
      smi2 = Chem.MolToSmiles(nms[0])
      self.assertEqual(csmi,smi2)
    ms = [Chem.MolFromSmiles(smi) for smi in smis]
    json = rdMolInterchange.MolsToJSON(ms)
    nms = rdMolInterchange.JSONToMols(json)
    self.assertEqual(len(ms),len(nms))
    self.assertEqual([Chem.MolToSmiles(x) for x in ms],[Chem.MolToSmiles(x) for x in nms])

  def test2(self):
    smis = ("C[C@H](O)C[C@@H](C)F |o1:1,4|","C[C@H](O)CC[C@@H](C)F |&1:1,5|")
    ms = [Chem.MolFromSmiles(x) for x in smis]
    json = rdMolInterchange.MolToJSON(ms[0])
    self.assertIn('stereoGroups',json)
    self.assertIn('"or"',json)
    self.assertIn('[1,4]',json)
    
    json = rdMolInterchange.MolToJSON(ms[1])
    self.assertIn('stereoGroups',json)
    self.assertIn('"and"',json)
    self.assertIn('[1,5]',json)
    
    json = rdMolInterchange.MolsToJSON(ms)
    self.assertIn('stereoGroups',json)
    self.assertIn('"or"',json)
    self.assertIn('[1,4]',json)
    self.assertIn('"and"',json)
    self.assertIn('[1,5]',json)
    
    ps = rdMolInterchange.JSONWriteParameters()
    ps.useRDKitExtensions = False
    json = rdMolInterchange.MolToJSON(ms[1],ps)
    self.assertNotIn('stereoGroups',json)
    json = rdMolInterchange.MolsToJSON(ms,ps)
    self.assertNotIn('stereoGroups',json)



  
if __name__ == '__main__':
  unittest.main()
