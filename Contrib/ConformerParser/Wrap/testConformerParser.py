from rdkit import Chem
from rdkit.Chem import rdConformerParser
from rdkit import RDConfig
import unittest
import os
from rdkit.RDLogger import logger
logger = logger()


class TestCase(unittest.TestCase):

  def setUp(self):
    pass

  def testReadAmberTraj(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Contrib', 'ConformerParser', 'test_data',
                         'water_coords.trx')
    mol = Chem.MolFromSmiles('O')
    mol = Chem.AddHs(mol)
    ids = rdConformerParser.AddConformersFromAmberTrajectory(mol, fileN)
    self.failUnless(mol.GetNumConformers() == 1)
    self.failUnless(len(ids) == 1)
    self.failUnless(ids[0] == 0)

    fileN = os.path.join(RDConfig.RDBaseDir, 'Contrib', 'ConformerParser', 'test_data',
                         'water_coords2.trx')
    ids = rdConformerParser.AddConformersFromAmberTrajectory(mol, fileN, clearConfs=True)
    self.failUnless(mol.GetNumConformers() == 2)
    ids = rdConformerParser.AddConformersFromAmberTrajectory(mol, fileN, clearConfs=False)
    self.failUnless(mol.GetNumConformers() == 4)
    ids = rdConformerParser.AddConformersFromAmberTrajectory(mol, fileN, numConfs=1,
                                                             clearConfs=True)
    self.failUnless(mol.GetNumConformers() == 1)


if __name__ == '__main__':
  unittest.main()
