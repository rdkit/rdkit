import os
import unittest

from rdkit import Chem, RDConfig
from rdkit.Chem import rdConformerParser
from rdkit.RDLogger import logger

logger = logger()


class TestCase(unittest.TestCase):

  def setUp(self):
    pass

  def testReadAmberTraj(self):
    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'test_data', 'water_coords.trx')
    mol = Chem.MolFromSmiles('O')
    mol = Chem.AddHs(mol)
    ids = rdConformerParser.AddConformersFromAmberTrajectory(mol, fileN)
    self.assertTrue(mol.GetNumConformers() == 1)
    self.assertTrue(len(ids) == 1)
    self.assertTrue(ids[0] == 0)

    fileN = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'test_data', 'water_coords2.trx')
    ids = rdConformerParser.AddConformersFromAmberTrajectory(mol, fileN, clearConfs=True)
    self.assertTrue(mol.GetNumConformers() == 2)
    ids = rdConformerParser.AddConformersFromAmberTrajectory(mol, fileN, clearConfs=False)
    self.assertTrue(mol.GetNumConformers() == 4)
    ids = rdConformerParser.AddConformersFromAmberTrajectory(mol, fileN, numConfs=1,
                                                             clearConfs=True)
    self.assertTrue(mol.GetNumConformers() == 1)


if __name__ == '__main__':
  unittest.main()
