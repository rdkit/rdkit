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


if __name__ == '__main__':
  unittest.main()
