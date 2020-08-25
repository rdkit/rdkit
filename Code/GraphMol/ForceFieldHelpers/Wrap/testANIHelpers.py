from rdkit import Chem
from rdkit.Chem import ChemicalForceFields, rdDistGeom
from rdkit import RDConfig
import unittest
import os
import numpy


class TestCase(unittest.TestCase):

  def setUp(self):
    self.dirName = os.path.join(RDConfig.RDBaseDir, 'Code', 'GraphMol', 'ForceFieldHelpers', 'UFF',
                                'test_data')

  def testANIBuilder(self):
    fName = os.path.join(self.dirName, 'CH4.mol')
    m = Chem.MolFromMolFile(fName, True, False)

    ff = ChemicalForceFields.ANIGetMoleculeForceField(m, "ANI-1ccx", 8)
    self.failUnless(ff)
    positions = ff.Positions()
    savedPos = list(positions)
    e1 = ff.CalcEnergy(savedPos)

    self.failUnless((e1 - (-40.0553)) < 0.05)
    r = ff.Minimize()
    self.failUnless(r == 0)
    e2 = ff.CalcEnergy()
    self.failUnless(e2 < e1)

    m1 = Chem.MolFromMolFile(fName, True, False)
    ff1 = ChemicalForceFields.ANIGetMoleculeForceField(m1, "ANI-1x", 8)
    self.failUnless(ff1)
    positions1 = ff1.Positions()
    savedPos1 = list(positions1)
    e1 = ff1.CalcEnergy(savedPos1)

    self.failUnless((e1 - (-40.0517)) < 0.05)
    r = ff1.Minimize()
    self.failUnless(r == 0)
    e2 = ff1.CalcEnergy()
    self.failUnless(e2 < e1)
  
  def ANIOptimizeMoleculeConfsTest(self): 
    fName = os.path.join(self.dirName, 'CH4.mol')
    m2 = Chem.MolFromMolFile(fName, True, False)
    ff2 = ChemicalForceFields.ANIOptimizeMoleculeConfs(m2)
    self.assertEqual(len(ff2), 1)
    self.assertEqual(ff2[0][0], 0)
    self.failUnless(ff2[0][1] < -40.0553)

  def testANIOptimizeMolecule(self):
    fName = os.path.join(self.dirName, 'CH4.mol')
    m2 = Chem.MolFromMolFile(fName, True, False)
    ff = ChemicalForceFields.ANIGetMoleculeForceField(m2, "ANI-1ccx", 8)
    self.failUnless(ff)
    e1 = ff.CalcEnergy()
    res = ChemicalForceFields.ANIOptimizeMolecule(m2)
    self.assertEqual(res, 0)
    ff1 = ChemicalForceFields.ANIGetMoleculeForceField(m2, "ANI-1ccx", 8)
    self.failUnless(ff1)
    e2 = ff1.CalcEnergy()
    self.failUnless(e2 < e1)


if __name__ == '__main__':
  unittest.main()
