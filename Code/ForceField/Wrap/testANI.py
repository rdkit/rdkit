from rdkit import RDConfig
import sys, os
from time import sleep
from multiprocessing import Process, Value
import unittest
from rdkit import Chem
from rdkit.Chem import ChemicalForceFields
from rdkit.Chem import rdMolTransforms

class TestCase(unittest.TestCase):
  def testANIForceField(self):
    self.dirName = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'ForceFieldHelpers', 'UFF',
                                'test_data')
    fName = os.path.join(self.dirName, 'CH4.mol')
    m = Chem.MolFromMolFile(fName, True, False)
    ff = ChemicalForceFields.ANIGetMoleculeForceField(m, "ANI-1ccx", 8)
    self.failUnless(ff)
    positions = ff.Positions()
    savedPos = list(positions)
    e1 = ff.CalcEnergy()
    ff.AddANIAtomContrib([0, 0, 0, 1, 0], 1, 3, 5, 0, 8, "ANI-1ccx")
    e2 = ff.CalcEnergy(savedPos)
    self.assertAlmostEqual(38.0375 - abs(e2 - e1), 0.0748, 3)

if __name__ == '__main__':
  unittest.main()